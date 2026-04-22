# ─────────────────────────────────────────────────
# src/Physics.jl
# ─────────────────────────────────────────────────

mass(p::Particle)   = SPECIES[p.species].mass
radius(p::Particle) = SPECIES[p.species].radius
kinetic_energy(p::Particle) = 0.5 * mass(p) * dot(p.vel, p.vel)
speed(p::Particle)          = norm(p.vel)

function wrap!(p::Particle, box::Box)
    p.pos[1] = mod(p.pos[1], box.width)
    p.pos[2] = mod(p.pos[2], box.height)
end

function min_image_delta(a::Particle, b::Particle, box::Box)
    dx = b.pos[1] - a.pos[1]
    dy = b.pos[2] - a.pos[2]
    dx -= box.width  * round(dx / box.width)
    dy -= box.height * round(dy / box.height)
    return dx, dy
end

"""
Two-phase hard-sphere collision resolver.

Phase 1 — positional correction:
  Move i and j apart along the line of centres until they just touch
  (overlap → 0).  Split in proportion to mass so the CM is conserved.
  This prevents the "sticky particle" bug where repeated impulses on the
  same overlapping pair accelerate them inward instead of outward.

Phase 2 — velocity impulse:
  Standard elastic collision formula in the CM frame, giving:
    Δv_i = −[2·m_j/(m_i+m_j)] (v_rel · n̂) n̂
    Δv_j = +[2·m_i/(m_i+m_j)] (v_rel · n̂) n̂
  Both total momentum and total kinetic energy are conserved exactly.
"""
function collide!(i::Particle, j::Particle, box::Box)
    dx, dy = min_image_delta(j, i, box)
    d2 = dx*dx + dy*dy

    sigma = radius(i) + radius(j)
    d2 > sigma^2  && return false
    d2 < 1e-14    && return false

    dvx = i.vel[1] - j.vel[1]
    dvy = i.vel[2] - j.vel[2]
    vdotn = dvx*dx + dvy*dy
    vdotn >= 0.0  && return false      # already separating

    # ── Phase 1: push apart ────────────────────────────────────────────────
    dist        = sqrt(d2)
    overlap     = sigma - dist
    inv_d       = 1.0 / dist
    nx, ny      = dx * inv_d, dy * inv_d
    mi, mj     = mass(i), mass(j)
    total_mass  = mi + mj
    fi          = mj / total_mass
    fj          = mi / total_mass
    i.pos[1]   += fi * overlap * nx;  i.pos[2] += fi * overlap * ny
    j.pos[1]   -= fj * overlap * nx;  j.pos[2] -= fj * overlap * ny
    wrap!(i, box);  wrap!(j, box)

    # ── Phase 2: velocity impulse ──────────────────────────────────────────
    factor      = 2.0 * mi * mj / total_mass * vdotn / d2
    i.vel[1]   -= factor / mi * dx;  i.vel[2] -= factor / mi * dy
    j.vel[1]   += factor / mj * dx;  j.vel[2] += factor / mj * dy

    return true
end

"""
Attempt the reaction A + A → B + C for two colliding A particles.

If the reaction fires (with probability `p_react`):

  Both total momentum AND total kinetic energy are conserved exactly.
  The derivation follows standard reactive-MD practice:

  1. Compute the pair's centre-of-mass (CM) velocity:
       v_cm = (p_i + p_j) / (m_A + m_A)
     This is invariant — the reaction cannot change it.

  2. Compute the kinetic energy available in the CM frame:
       KE_cm = ½ m_A |v_i - v_cm|² + ½ m_A |v_j - v_cm|²
             = ¼ m_A |v_rel|²

  3. Distribute KE_cm between products B and C via their reduced mass
       μ_BC = m_B m_C / (m_B + m_C):
       |v_rel'| = sqrt(2 KE_cm / μ_BC)
     The direction of v_rel' is chosen uniformly at random (isotropic
     scattering, no preferred axis after the reaction).

  4. Reconstruct lab-frame velocities:
       v_B = v_cm + [m_C/(m_B+m_C)] v_rel'
       v_C = v_cm − [m_B/(m_B+m_C)] v_rel'

  By construction: m_B v_B + m_C v_C = (m_B+m_C) v_cm = (2 m_A) v_cm = p_i + p_j  ✓
                   ½ m_B |v_B-v_cm|² + ½ m_C |v_C-v_cm|² = KE_cm  ✓

Returns true if the reaction fired.
"""
function react_AA!(i::Particle, j::Particle, p_react::Float64)
    # Solo aplica a pares A + A
    (i.species == 1 && j.species == 1) || return false
    rand() > p_react && return false

    m_A = SPECIES[1].mass
    m_B = SPECIES[2].mass
    m_C = SPECIES[3].mass

    M_old = 2.0 * m_A
    M_new = m_B + m_C

    # ── 1. Conservación estricta del Momento Total ─────────────────────────
    # Calculamos el momento absoluto previo a la reacción
    P_x = m_A * i.vel[1] + m_A * j.vel[1]
    P_y = m_A * i.vel[2] + m_A * j.vel[2]

    # La nueva velocidad del CM debe ajustarse si la masa cambió
    vcm_x_new = P_x / M_new
    vcm_y_new = P_y / M_new

    # ── 2. Conservación estricta de la Energía Cinética ────────────────────
    KE_total = 0.5 * m_A * (i.vel[1]^2 + i.vel[2]^2) + 
               0.5 * m_A * (j.vel[1]^2 + j.vel[2]^2)

    # Energía que "consume" el desplazamiento del nuevo centro de masas
    KE_cm_new = 0.5 * M_new * (vcm_x_new^2 + vcm_y_new^2)
    
    # Energía restante disponible para la repulsión de los productos
    KE_rel = KE_total - KE_cm_new

    # Guardia termodinámica: si la masa aumenta tanto que requiere más energía 
    # de la que existe, la reacción es "muy endotérmica" y no puede ocurrir.
    KE_rel < 0.0 && return false

    # ── 3. Reconstruir velocidades de los productos B y C ──────────────────
    mu_BC   = m_B * m_C / M_new
    vrel_cm = sqrt(2.0 * KE_rel / mu_BC)

    theta = 2π * rand()
    rx = vrel_cm * cos(theta)
    ry = vrel_cm * sin(theta)

    fB = m_C / M_new
    fC = m_B / M_new

    i.vel[1] = vcm_x_new + fB * rx;  i.vel[2] = vcm_y_new + fB * ry
    j.vel[1] = vcm_x_new - fC * rx;  j.vel[2] = vcm_y_new - fC * ry

    # ── 4. Cambiar especies ────────────────────────────────────────────────
    i.species = 2   # becomes B
    j.species = 3   # becomes C

    return true
end

function react_acid_base!(i::Particle, j::particle, p_react::Float64)
    # comprueba que sean el tipo de particulas correcta
    if ((i.species == 2 && j.species == 3) || (i.species == 3 && j.species == 2)) || return false
    
    if (i.species == 3 && j.species == 2)
        k = i
        i = j
        j = k
    end

    # probabilidad de reaccionar
    rand() > p_react && return false

    #-------------------- actualizar valores---------------------------
    m_h2o = SPECIES[1].mass
    m_h3o = SPECIES[2].mass
    m_oh = SPECIES[3].mass

    M_old = m_h3o + m_oh
    M_new = 2*m_h2o

    # ── 1. Conservación estricta del Momento Total ─────────────────────────
    # Calculamos el momento absoluto previo a la reacción
    P_x = m_h30 * i.vel[1] + m_oh * j.vel[1]
    P_y = m_h30 * i.vel[2] + m_oh * j.vel[2]

    # La nueva velocidad del CM debe ajustarse si la masa cambió
    vcm_x_new = P_x / M_new
    vcm_y_new = P_y / M_new

    # ── 2. Conservación estricta de la Energía Cinética ────────────────────
    KE_total = 0.5 * m_h30 * (i.vel[1]^2 + i.vel[2]^2) + 
               0.5 * m_oh * (j.vel[1]^2 + j.vel[2]^2)

    # Energía que "consume" el desplazamiento del nuevo centro de masas
    KE_cm_new = 0.5 * M_new * (vcm_x_new^2 + vcm_y_new^2)
    
    # Energía restante disponible para la repulsión de los productos
    KE_rel = KE_total - KE_cm_new

    # Guardia termodinámica: si la masa aumenta tanto que requiere más energía 
    # de la que existe, la reacción es "muy endotérmica" y no puede ocurrir.
    KE_rel < 0.0 && return false

    # ── 3. Reconstruir velocidades de los productos B y C ──────────────────
    mu_BC   = m_B * m_C / M_new
    vrel_cm = sqrt(2.0 * KE_rel / mu_BC)

    theta = 2π * rand()
    rx = vrel_cm * cos(theta)
    ry = vrel_cm * sin(theta)

    fB = m_C / M_new
    fC = m_B / M_new

    i.vel[1] = vcm_x_new + fB * rx;  i.vel[2] = vcm_y_new + fB * ry
    j.vel[1] = vcm_x_new - fC * rx;  j.vel[2] = vcm_y_new - fC * ry

    # ── 4. Cambiar especies ────────────────────────────────────────────────
    i.species = 2   # becomes B
    j.species = 3   # becomes C

end

"""
One simulation step:
  1. Detect and resolve all pairwise overlaps.
     For A+A collisions, attempt the chemical reaction first;
     if it does not fire, apply the elastic collision.
  2. Advance positions with velocity-Verlet (no forces, so just free streaming).
  3. Apply periodic boundary conditions.

Returns (n_elastic, n_reactions).
"""
function step!(particles::Vector{Particle}, box::Box, dt::Float64, p_react::Float64)
    n          = length(particles)
    n_elastic  = 0
    n_reactions = 0

    for ii in 1:n-1, jj in ii+1:n
        pi = particles[ii]
        pj = particles[jj]

        # Check geometric overlap first (cheap)
        dx, dy = min_image_delta(pj, pi, box)
        d2 = dx*dx + dy*dy
        sigma = radius(pi) + radius(pj)
        d2 > sigma^2 && continue
        d2 < 1e-14   && continue

        # A+A pair: try reaction before elastic bounce
        if pi.species == 1 && pj.species == 1
            if react_AA!(pi, pj, p_react)
                n_reactions += 1
                # No positional correction — see collide! docstring.
                continue
            end
        end

        # Standard elastic collision
        if collide!(pi, pj, box)
            n_elastic += 1
        end
    end

    # Free streaming
    for p in particles
        p.pos[1] += p.vel[1] * dt
        p.pos[2] += p.vel[2] * dt
        wrap!(p, box)
    end

    return n_elastic, n_reactions
end

function diagnostics(particles::Vector{Particle})
    KE   = sum(kinetic_energy(p) for p in particles)
    px   = sum(mass(p) * p.vel[1] for p in particles)
    py   = sum(mass(p) * p.vel[2] for p in particles)
    T    = KE / length(particles)
    N_A  = count(p -> p.species == 1, particles)
    N_B  = count(p -> p.species == 2, particles)
    N_C  = count(p -> p.species == 3, particles)
    return (; KE, px, py, T, N_A, N_B, N_C)
end