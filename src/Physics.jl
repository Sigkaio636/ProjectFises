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
Attempt a stochastic binary chemical reaction between two particles, mutating
them in-place if the reaction succeeds. Returns `true` if the reaction occurred,
`false` otherwise.

The reaction modelled is of the form:

    Aᵢ + Aⱼ  ──p_react──►  Bᵢ + Bⱼ

where the input species `(i.species, j.species)` are transmuted into the output
species `(iout_sp, jout_sp)`. The function enforces **exact conservation** of
linear momentum and kinetic energy (with a chemical energy correction) at every
accepted event.

# Arguments
- `i::Particle`       : First reactant particle (mutated in-place on success).
- `j::Particle`       : Second reactant particle (mutated in-place on success).
- `iout_sp::Int`      : Species index (into `SPECIES`) assigned to `i` after the reaction.
- `jout_sp::Int`      : Species index (into `SPECIES`) assigned to `j` after the reaction.
- `p_react::Float64`  : Bare stochastic probability of reaction per collision
                        attempt, ∈ [0, 1]. Acts as an effective reaction rate
                        constant at the collision-detection level.

# Returns
- `true`  — reaction occurred; `i` and `j` have been mutated.
- `false` — reaction did not occur (stochastic rejection **or** thermodynamic
            veto); `i` and `j` are unchanged.

# Physics & Algorithm

## 1 — Stochastic gate
A uniform random number is drawn first. If it exceeds `p_react` the function
returns immediately with no side-effects, modelling a reaction probability per
collision.

## 2 — Total momentum conservation
The total linear momentum **P** = mᵢ **vᵢ** + mⱼ **vⱼ** is computed before the
reaction. Because the total mass may change (Mₒₗd ≠ Mₙₑw), the centre-of-mass
(CM) velocity is recomputed from the *new* total mass so that **P** is
preserved exactly:

    vcm_new = P / M_new

## 3 — Energy balance
Total kinetic energy before the reaction is split into:

    E_init = KE_cm_new + KE_rel + E_react

| Term        | Meaning                                                      |
|:----------- |:-------------------------------------------------------------|
| `E_init`    | Total pre-reaction kinetic energy of the two particles       |
| `KE_cm_new` | KE of the new CM motion (fixed by momentum conservation)     |
| `E_react`   | Internal chemical energy change  (ΔE = Σeₒᵤₜ − Σeᵢₙ)       |
| `KE_rel`    | Remaining energy available for relative motion of products   |

`E_react > 0` is endothermic (energy absorbed); `E_react < 0` is exothermic
(energy released into relative motion).

## 4 — Thermodynamic veto
If `KE_rel < 0` the reaction is forbidden: the collision does not carry enough
kinetic energy to supply the endothermic demand (or the mass increase), and the
function returns `false` without modifying the particles. This is the analogue
of an activation-energy threshold.

## 5 — Product velocity reconstruction
The relative speed of the products is obtained from the reduced mass μ = mᵢₒᵤₜ mⱼₒᵤₜ / Mₙₑw:

    v_rel = √(2 KE_rel / μ)

The direction θ is drawn uniformly on [0, 2π), modelling an isotropic
(hard-sphere-like) collision in 2-D. Velocities are then reconstructed in the
lab frame:

    vᵢ = vcm_new + (mⱼₒᵤₜ / Mₙₑw) * v_rel_vec
    vⱼ = vcm_new − (mᵢₒᵤₜ / Mₙₑw) * v_rel_vec

## 6 — Species mutation
Only after all checks pass are `i.species` and `j.species` updated to
`iout_sp` and `jout_sp` respectively.

# Conservation laws verified
| Quantity              | Conserved? | Notes                                      |
|:----------------------|:----------:|:-------------------------------------------|
| Linear momentum (x,y) | ✅ exact   | holds even when M_old ≠ M_new              |
| Kinetic + chemical E  | ✅ exact   | reaction vetoed if budget is insufficient  |
| Angular momentum      | ❌         | not tracked (point-particle 2-D model)     |
| Particle count        | ✅         | two in, two out                            |

# Side-effects
Mutates `i.vel`, `j.vel`, `i.species`, and `j.species` **only** on `return true`.
The global `SPECIES` table is read but never modified.

# Errors / assumptions
- `SPECIES` must be accessible in scope and indexed by the integer species codes.
- Species indices `iout_sp` and `jout_sp` must be valid keys in `SPECIES`.
- `p_react` is not validated; passing values outside [0,1] produces undefined
  stochastic behaviour.
- The function is **not** thread-safe if `i` and `j` can be accessed concurrently.
"""

function react_chem!(i::Particle, j::Particle, iout_sp::Int, jout_sp::Int, p_react::Float64)
    rand() > p_react && return false

    m_i = SPECIES[i.species].mass
    m_j = SPECIES[j.species].mass
    m_iout = SPECIES[iout_sp].mass
    m_jout = SPECIES[jout_sp].mass

    e_i = SPECIES[i.species].energy
    e_j = SPECIES[j.species].energy
    e_iout = SPECIES[iout_sp].energy
    e_jout = SPECIES[jout_sp].energy

    M_old = m_i + m_j
    M_new = m_iout + m_jout

    # ── 1. Conservación estricta del Momento Total ─────────────────────────
    # Calculamos el momento absoluto previo a la reacción
    P_x = m_i * i.vel[1] + m_j * j.vel[1]
    P_y = m_i * i.vel[2] + m_j * j.vel[2]

    # La nueva velocidad del CM debe ajustarse si la masa cambió
    vcm_x_new = P_x / M_new
    vcm_y_new = P_y / M_new

    # ── 2. Conservación estricta de la Energía ────────────────────
    # Al principio es solo energía cinética de las particulas incidentes
    E_init = 0.5 * m_i * (i.vel[1]^2 + i.vel[2]^2) + 
               0.5 * m_j * (j.vel[1]^2 + j.vel[2]^2)

    # Energía que "consume" el desplazamiento del nuevo centro de masas
    KE_cm_new = 0.5 * M_new * (vcm_x_new^2 + vcm_y_new^2)

    # Energía de reacción entre las particulas
    E_react  = e_iout + e_jout - e_i - e_j

    # Balance de energía: E_init = KE_rel + KE_cm_new + E_react
    
    # Energía restante disponible para la repulsión de los productos
    KE_rel = E_init - KE_cm_new - E_react

    # Guardia termodinámica: si la reacción es "muy endotérmica" 
    #   (porque masa aumenta tanto que requiere más energía de la que existe,
    #    o porque la reacción química posee una gran energía de activación)
    #   y no puede ocurrir.
    KE_rel < 0.0 && return false

    # ── 3. Reconstruir velocidades de los productos B y C ──────────────────
    mu_BC   = m_iout * m_jout / M_new
    vrel_cm = sqrt(2.0 * KE_rel / mu_BC)

    theta = 2π * rand()
    rx = vrel_cm * cos(theta)
    ry = vrel_cm * sin(theta)

    fiout = m_jout / M_new
    fjout = m_iout / M_new

    i.vel[1] = vcm_x_new + fiout * rx;  i.vel[2] = vcm_y_new + fiout * ry
    j.vel[1] = vcm_x_new - fjout * rx;  j.vel[2] = vcm_y_new - fjout * ry

    # ── 4. Cambiar especies ────────────────────────────────────────────────
    i.species = iout_sp
    j.species = jout_sp

    return true
end

"""
One simulation step:
  1. Detect and resolve all pairwise overlaps.
     For A+A collisions, attempt the chemical reaction first;
     if it does not fire, apply the elastic collision.
  2. Advance positions with velocity-Verlet (no forces, so just free streaming).
  3. Apply periodic boundary conditions.

Returns (n_collitions, n_reactions).
"""
function step!(particles::Vector{Particle}, box::Box, dt::Float64, p_react::Float64)
    n          = length(particles)
    n_collitions  = 0
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

        # H2O+H2O pair: try reaction before elastic bounce
        if pi.species == 1 && pj.species == 1
            if react_chem!(pi, pj, 2, 3, p_react)
                n_reactions += 1
                # No positional correction — see collide! docstring.
                continue
            end
        end

        # Standard elastic collision
        if collide!(pi, pj, box)
            n_collitions += 1
        end
    end

    # Free streaming
    for p in particles
        p.pos[1] += p.vel[1] * dt
        p.pos[2] += p.vel[2] * dt
        wrap!(p, box)
    end

    return n_collitions, n_reactions
end

function diagnostics(particles::Vector{Particle})
    KE   = sum(kinetic_energy(p) for p in particles)
    px   = sum(mass(p) * p.vel[1] for p in particles)
    py   = sum(mass(p) * p.vel[2] for p in particles)
    T    = 2/2*KE / length(particles)
    N_H2O  = count(p -> p.species == 1, particles)
    N_H3O  = count(p -> p.species == 2, particles)
    N_OH  = count(p -> p.species == 3, particles)
    return (; KE, px, py, T, N_H2O, N_H3O, N_OH)
end