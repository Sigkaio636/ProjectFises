# ─────────────────────────────────────────────────
# src/Physics.jl
# ─────────────────────────────────────────────────

kinetic_energy(p::Particle) = 0.5 * p.mass * dot(p.vel, p.vel)
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

function collide!(i::Particle, j::Particle, box::Box)
    dx, dy = min_image_delta(j, i, box)
    d2 = dx*dx + dy*dy

    sigma = i.radius + j.radius
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
    total_mass  = i.mass + j.mass
    fi          = j.mass / total_mass
    fj          = i.mass / total_mass
    i.pos[1]   += fi * overlap * nx;  i.pos[2] += fi * overlap * ny
    j.pos[1]   -= fj * overlap * nx;  j.pos[2] -= fj * overlap * ny
    wrap!(i, box);  wrap!(j, box)

    # ── Phase 2: velocity impulse ──────────────────────────────────────────
    factor      = 2.0 * i.mass * j.mass / total_mass * vdotn / d2
    i.vel[1]   -= factor / i.mass * dx;  i.vel[2] -= factor / i.mass * dy
    j.vel[1]   += factor / j.mass * dx;  j.vel[2] += factor / j.mass * dy

    return true
end

function step!(particles::Vector{Particle}, box::Box, dt::Float64)
    n = length(particles)
    ncollisions = 0
    for i in 1:n-1, j in i+1:n
        collide!(particles[i], particles[j], box) && (ncollisions += 1)
    end
    for p in particles
        p.pos[1] += p.vel[1] * dt
        p.pos[2] += p.vel[2] * dt
        wrap!(p, box)
    end
    return ncollisions
end

function diagnostics(particles::Vector{Particle})
    KE  = sum(kinetic_energy(p) for p in particles)
    px  = sum(p.mass * p.vel[1] for p in particles)
    py  = sum(p.mass * p.vel[2] for p in particles)
    T   = KE / length(particles)
    return (; KE, px, py, T)
end