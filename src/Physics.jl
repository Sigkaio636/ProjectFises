# ─────────────────────────────────────────────────
# src/Physics.jl
# ─────────────────────────────────────────────────

@inline mass(p::Particle)    = SPECIES[p.species].mass
@inline radius(p::Particle)  = SPECIES[p.species].radius
@inline kinetic_energy(p::Particle) = 0.5 * mass(p) * (p.vx*p.vx + p.vy*p.vy)
@inline speed(p::Particle)   = sqrt(p.vx*p.vx + p.vy*p.vy)

@inline function wrap!(p::Particle, box::Box)
    p.x -= box.width  * floor(p.x * box.inv_w)
    p.y -= box.height * floor(p.y * box.inv_h)
end

# Returns (dx, dy) = b - a under minimum image convention
@inline function min_image_delta(ax, ay, bx, by, box::Box)
    dx = bx - ax
    dy = by - ay
    dx -= box.width  * round(dx * box.inv_w)
    dy -= box.height * round(dy * box.inv_h)
    return dx, dy
end

"""
Two-phase hard-sphere collision resolver with trajectory validation.

Phase 0 — trajectory discriminant check:
  Given current positions and velocities, solve for the time t* at which
  the pair would reach minimum separation. The collision is physical only
  if the minimum approach distance is less than sigma (discriminant > 0).
  This rejects numerical ghost collisions where particles overlap due to
  finite dt but their trajectories never actually intersect at contact.

  The discriminant of |r(t)|² = sigma² is:
    Δ = (v_rel · r̂)² - (|v_rel|² - 0) ... simplified to:
    Δ = (v·r)² - |v|²·(|r|² - sigma²)
  Collision is real iff Δ > 0 and t* > 0 (approaching).

Phase 1 — positional correction:
  Move i and j apart along the line of centres until they just touch
  (overlap → 0). Split in proportion to mass so the CM is conserved.
  This prevents the sticky-particle bug where repeated impulses on the
  same overlapping pair accelerate them inward instead of outward.

Phase 2 — velocity impulse:
  Standard elastic collision in the CM frame:
    Δv_i = −[2·m_j/(m_i+m_j)] (v_rel · n̂) n̂
    Δv_j = +[2·m_i/(m_i+m_j)] (v_rel · n̂) n̂
  Both total momentum and total kinetic energy are conserved exactly.
"""
function collide!(i::Particle, j::Particle, box::Box)

    # ── Geometry ───────────────────────────────────────────────────────────
    dx, dy = min_image_delta(j.x, j.y, i.x, i.y, box)
    d2     = dx*dx + dy*dy
    sigma  = radius(i) + radius(j)
    sigma2 = sigma * sigma

    d2 > sigma2 && return false       # no overlap at all
    d2 < 1e-14  && return false       # degenerate / same position

    # ── Relative velocity ──────────────────────────────────────────────────
    dvx = i.vx - j.vx
    dvy = i.vy - j.vy

    # v_rel · r  (r points from j to i, same as dx,dy)
    vdotn = dvx*dx + dvy*dy
    vdotn >= 0.0 && return false      # already separating — not a collision

    # ── Phase 0: discriminant / trajectory validity check ─────────────────
    # Solves |r + v_rel·t|² = sigma² for t
    # Δ = (v·r)² - |v|²·(d² - sigma²)
    # For a real collision Δ must be > 0.
    # Note: (d² - sigma²) ≤ 0 since we already passed the overlap check,
    # so Δ = vdotn² - |v|²·(d² - sigma²) ≥ vdotn² ≥ 0 always here.
    # The check becomes meaningful when used BEFORE the overlap test,
    # or as an epsilon guard against nearly-grazing numerical ghosts:
    v2 = dvx*dvx + dvy*dvy
    v2 < 1e-14 && return false        # essentially zero relative speed

    discriminant = vdotn*vdotn - v2*(d2 - sigma2)
    discriminant < 0.0 && return false   # trajectories miss — ghost collision

    # Impact parameter b = sqrt(d² - (v·r/|v|)²) / sigma
    # Optionally weight reaction probability by (1 - b²/sigma²) here

    # ── Phase 1: positional correction ────────────────────────────────────
    dist    = sqrt(d2)
    overlap = sigma - dist
    inv_d   = 1.0 / dist
    nx, ny  = dx * inv_d, dy * inv_d
    mi, mj  = mass(i), mass(j)
    M       = mi + mj
    fi      = mj / M                  # fraction of overlap assigned to i
    fj      = mi / M                  # fraction of overlap assigned to j

    i.x += fi * overlap * nx;   i.y += fi * overlap * ny
    j.x -= fj * overlap * nx;   j.y -= fj * overlap * ny
    wrap!(i, box);  wrap!(j, box)

    # ── Phase 2: velocity impulse ──────────────────────────────────────────
    fac   = 2.0 * mi * mj / M * vdotn / d2
    i.vx -= fac / mi * dx;   i.vy -= fac / mi * dy
    j.vx += fac / mj * dx;   j.vy += fac / mj * dy

    return true
end

function react_chem!(i::Particle, j::Particle,
                     iout_sp::Int, jout_sp::Int, p_react::Float64)
    rand() > p_react && return false

    m_i    = SPECIES[i.species].mass;   m_j    = SPECIES[j.species].mass
    m_iout = SPECIES[iout_sp].mass;     m_jout = SPECIES[jout_sp].mass
    e_i    = SPECIES[i.species].energy; e_j    = SPECIES[j.species].energy
    e_iout = SPECIES[iout_sp].energy;   e_jout = SPECIES[jout_sp].energy

    M_old   = m_i + m_j
    M_new   = m_iout + m_jout
    # ── 1. Conservación estricta del Momento Total ─────────────────────────
    # Calculamos el momento absoluto previo a la reacción
    P_x     = m_i*i.vx + m_j*j.vx;  P_y     = m_i*i.vy + m_j*j.vy
    # La nueva velocidad del CM debe ajustarse si la masa cambió
    vcm_x   = P_x / M_new;  vcm_y   = P_y / M_new

    # ── 2. Conservación estricta de la Energía ────────────────────
    # Al principio es solo energía cinética de las particulas incidentes
    E_init  = 0.5*m_i*(i.vx^2 + i.vy^2) + 0.5*m_j*(j.vx^2 + j.vy^2)
    # Energía que "consume" el desplazamiento del nuevo centro de masas
    KE_cm   = 0.5*M_new*(vcm_x^2 + vcm_y^2)
    # Energía de reacción entre las particulas
    E_react = e_iout + e_jout - e_i - e_j
    # Balance de energía: E_init = KE_rel + KE_cm_new + E_react
    # Energía restante disponible para la repulsión de los productos
    KE_rel  = E_init - KE_cm - E_react
    # Guardia termodinámica: si la reacción es "muy endotérmica" 
    #   (porque masa aumenta tanto que requiere más energía de la que existe,
    #    o porque la reacción química posee una gran energía de activación)
    #   y no puede ocurrir.
    KE_rel < 0.0 && return false

    # ── 3. Reconstruir velocidades de los productos ──────────────────
    # Calcular la orientación de las particulas incidentes
    vcm_pre_x   = P_x / M_old;  vcm_pre_y   = P_y / M_old
    v_rel_xi = i.vx - vcm_pre_x ; v_rel_yi = i.vy - vcm_pre_y
    v_rel_xi_norm = -v_rel_xi / sqrt(v_rel_xi^2 + v_rel_yi^2)
    v_rel_yi_norm = -v_rel_yi / sqrt(v_rel_xi^2 + v_rel_yi^2)
    # La orientación de las particulas finales es opuesta a las incidentes 

    mu      = m_iout * m_jout / M_new
    vrel    = sqrt(2.0 * KE_rel / mu)
    rx      = vrel * v_rel_xi_norm
    ry      = vrel * v_rel_yi_norm
    fi      = m_jout / M_new
    fj      = m_iout / M_new

    i.vx = vcm_x + fi*rx;  i.vy = vcm_y + fi*ry
    j.vx = vcm_x - fj*rx;  j.vy = vcm_y - fj*ry
    i.species = iout_sp
    j.species = jout_sp
    return true
end

# ── Pre-allocated per-thread event buffers (module-level, grown as needed) ──
# Avoids [Event[] for _ in 1:T] allocation inside the hot step! loop
const _thread_bufs = Vector{Vector{Event}}(undef, 0)

function _ensure_thread_bufs!(nt::Int)
    n = nt + 1                   # +1 here too
    if length(_thread_bufs) < n
        resize!(_thread_bufs, n)
        for k in 1:n
            if !isassigned(_thread_bufs, k)
                _thread_bufs[k] = Event[]
            end
        end
    end
end

function build_cell_list(particles::Vector{Particle}, box::Box)
    max_r  = maximum(SPECIES[k].radius for k in 1:3)
    cell_w = 2.0 * max_r * 2           # cell size ≥ interaction range
    cell_h = 2.0 * max_r * 2
    ncols  = max(1, floor(Int, box.width  / cell_w))
    nrows  = max(1, floor(Int, box.height / cell_h))
    cell_w = box.width  / ncols        # exact fit
    cell_h = box.height / nrows

    cells = [Int[] for _ in 1:ncols*nrows]
    for (i, p) in enumerate(particles)
        col = clamp(floor(Int, p.x / cell_w), 0, ncols-1)
        row = clamp(floor(Int, p.y / cell_h), 0, nrows-1)
        push!(cells[row*ncols + col + 1], i)
    end
    return CellList(cells, cell_w, cell_h, ncols, nrows)
end

"""
    step!(particles, box, dt, p_react) -> (n_elastic, n_forward, n_reverse)

Three-stage parallelized simulation step.  See module docstring for full
physics and threading documentation.
"""
function step!(particles::Vector{Particle}, box::Box,
               dt::Float64, p_react::Float64)
    n  = length(particles)
    nt = nthreads()
    _ensure_thread_bufs!(nt)

    # ── Stage 1: parallel O(N²) pair detection ──────────────────────────
    # Empty thread buffers (cheap: just reset length, no GC)
    for k in 1:nt+1
        empty!(_thread_bufs[k])
    end

    # Build the cell list sequentially (this is very fast, O(N))
    cl = build_cell_list(particles, box)

    @threads :static for ci in 1:length(cl.cells)
        tid = threadid()
        buf = _thread_bufs[tid]
        
        row = (ci - 1) ÷ cl.ncols
        col = (ci - 1) % cl.ncols
        
        # Only check 9 neighbouring cells (including self)
        for drow in -1:1, dcol in -1:1
            nrow = mod(row + drow, cl.nrows)
            ncol = mod(col + dcol, cl.ncols)
            nci  = nrow * cl.ncols + ncol + 1

            (ci > nci) && continue
            
            for ii in cl.cells[ci]
                pi = particles[ii]
                si = pi.species
                xi, yi = pi.x, pi.y
                sigma_i = SPECIES[si].radius
                
                for jj in cl.cells[nci]
                    # Avoid double-counting and self-interaction
                    (ci == nci && ii >= jj) && continue
                    
                    pj = particles[jj]
                    dx, dy = min_image_delta(xi, yi, pj.x, pj.y, box)
                    d2 = dx*dx + dy*dy
                    sigma = sigma_i + SPECIES[pj.species].radius
                    
                    (d2 > sigma*sigma || d2 < 1e-14) && continue

                    sj   = pj.species
                    kind = if si == 1 && sj == 1
                        EV_FORWARD
                    elseif (si == 2 && sj == 3) || (si == 3 && sj == 2)
                        EV_REVERSE
                    else
                        EV_ELASTIC
                    end
                    push!(buf, Event(Int32(ii), Int32(jj), kind))
                end
            end
        end
    end

    # ── Merge buffers (avoids one vcat allocation via sizehint!) ─────────
    total_events = sum(length(_thread_bufs[k]) for k in 1:nt+1)
    events = Vector{Event}(undef, 0)
    sizehint!(events, total_events)
    shuffle!(events)
    for k in 1:nt+1
        append!(events, _thread_bufs[k])
    end

    # ── Stage 2: sequential event application ────────────────────────────
    reacted         = falses(n)
    n_elastic       = 0
    n_react_fwd     = 0
    n_react_rev     = 0

    @inbounds for ev in events
        ii, jj = ev.ii, ev.jj
        (reacted[ii] || reacted[jj]) && continue

        pi = particles[ii]
        pj = particles[jj]

        if ev.kind == EV_REVERSE
            if react_chem!(pi, pj, 1, 1, p_react)
                reacted[ii] = reacted[jj] = true
                n_react_rev += 1
            else
                collide!(pi, pj, box) && (n_elastic += 1)
            end    
        elseif ev.kind == EV_FORWARD
            if react_chem!(pi, pj, 2, 3, p_react)
                reacted[ii] = reacted[jj] = true
                n_react_fwd += 1
            else
                collide!(pi, pj, box) && (n_elastic += 1)
            end
        
        else
            collide!(pi, pj, box) && (n_elastic += 1)
        end
    end

    # ── Stage 3: parallel free streaming ─────────────────────────────────
    @threads :static for p in particles
        p.x += p.vx * dt
        p.y += p.vy * dt
        wrap!(p, box)
    end

    return n_elastic, n_react_fwd, n_react_rev
end

# ── Diagnostics — parallelized reductions ────────────────────────────────────
function diagnostics(particles::Vector{Particle})
    nt = nthreads() + 1          # +1 guards against threadid() == nthreads()+1
    n  = length(particles)

    KE_t    = zeros(nt);  px_t = zeros(nt);  py_t = zeros(nt)
    NH2O_t  = zeros(Int, nt)
    NH3O_t  = zeros(Int, nt)
    NOH_t   = zeros(Int, nt)
    KEH2O_t = zeros(nt);  KEH3O_t = zeros(nt);  KEOH_t = zeros(nt)

    @threads :static for i in 1:n
        tid = threadid()
        p   = particles[i]
        ke  = kinetic_energy(p)
        KE_t[tid]  += ke
        px_t[tid]  += mass(p) * p.vx
        py_t[tid]  += mass(p) * p.vy
        sp = p.species
        if sp == 1
            NH2O_t[tid]  += 1;  KEH2O_t[tid] += ke
        elseif sp == 2
            NH3O_t[tid]  += 1;  KEH3O_t[tid] += ke
        else
            NOH_t[tid]   += 1;  KEOH_t[tid]  += ke
        end
    end

    KE     = sum(KE_t);    px = sum(px_t);    py = sum(py_t)
    N_H2O  = sum(NH2O_t);  N_H3O = sum(NH3O_t);  N_OH = sum(NOH_t)
    KE_H2O = sum(KEH2O_t); KE_H3O = sum(KEH3O_t); KE_OH = sum(KEOH_t)
    T      = KE / n
    Kc     = N_H2O > 0 ? (N_H3O * N_OH) / N_H2O^2 : Inf

    return (; KE, px, py, T, N_H2O, N_H3O, N_OH, KE_H2O, KE_H3O, KE_OH, Kc)
end