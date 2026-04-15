# ─────────────────────────────────────────────────
# src/Types.jl
# ─────────────────────────────────────────────────

mutable struct Particle
    pos::Vector{Float64}   # [x, y]
    vel::Vector{Float64}   # [vx, vy]
    mass::Float64
    radius::Float64
end

Particle(pos, vel; mass=1.0, radius=0.5) = Particle(pos, vel, mass, radius)

struct Box
    width::Float64
    height::Float64
end