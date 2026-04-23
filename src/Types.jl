# ─────────────────────────────────────────────────
# src/Types.jl
# ─────────────────────────────────────────────────

const SPECIES = (
    A = (mass=1.0, radius=2.0, energy=5, color=:mediumseagreen, name="H2O"),
    B = (mass=1.0, radius=2.0, energy=0, color=:tomato, name="H3O+"),
    C = (mass=1.0, radius=2.0, energy=0, color=:cornflowerblue, name="OH-"),
)

mutable struct Particle
    pos::Vector{Float64}   # [x, y]
    vel::Vector{Float64}   # [vx, vy]
    species::Int           # 1=H2O  2=H3O+  3=OH-
end

function Particle(pos, vel, species::Int)
    @assert species in (1, 2, 3) "species must be 1, 2, or 3"
    Particle(pos, vel, species)
end

struct Box
    width::Float64
    height::Float64
end