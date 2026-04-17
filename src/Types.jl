# ─────────────────────────────────────────────────
# src/Types.jl
# ─────────────────────────────────────────────────

"""
Each species has a fixed mass, radius, and display colour.
Species indices: 1 = A, 2 = B, 3 = C
"""
const SPECIES = (
    A = (mass=1.0, radius=1.0, color=:tomato,         name="A"),
    B = (mass=2.5, radius=2.8, color=:mediumseagreen, name="B"), # 1.5 + 0.5 = 2.0
    C = (mass=3.5, radius=3.4, color=:cornflowerblue, name="C"),
)

mutable struct Particle
    pos::Vector{Float64}   # [x, y]
    vel::Vector{Float64}   # [vx, vy]
    species::Int           # 1=A  2=B  3=C
end

function Particle(pos, vel, species::Int)
    @assert species in (1, 2, 3) "species must be 1, 2, or 3"
    Particle(pos, vel, species)
end

struct Box
    width::Float64
    height::Float64
end