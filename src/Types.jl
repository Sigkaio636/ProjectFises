# ─────────────────────────────────────────────────
# src/Types.jl
# ─────────────────────────────────────────────────

# Immutable species table — accessed read-only from all threads
struct SpeciesData
    mass    :: Float64
    radius  :: Float64
    energy  :: Float64
    color   :: Symbol
    name    :: String
end

const SPECIES = (
    SpeciesData(1.0, 2.0, 2.0, :mediumseagreen,   "H₂O" ),   # 1
    SpeciesData(1.0, 2.0, 0.0, :tomato,            "H₃O⁺"),   # 2
    SpeciesData(1.0, 2.0, 0.0, :cornflowerblue,    "OH⁻" ),   # 3
)

# SoA-friendly flat struct: all fields are isbits → stack-allocated, no GC pressure
mutable struct Particle
    x  :: Float64
    y  :: Float64
    vx :: Float64
    vy :: Float64
    species :: Int32          # Int32 saves 4 bytes, fits in same cache line
end

struct Box
    width  :: Float64
    height :: Float64
    inv_w  :: Float64         # precomputed reciprocals avoid division in wrap!
    inv_h  :: Float64
end
Box(w, h) = Box(w, h, 1.0/w, 1.0/h)

# Encode kind as Int8 to keep Event to 3×8 + 1 = 17 → padded to 24 bytes
# Using constants avoids symbol hash lookups in hot loop
const EV_ELASTIC  = Int8(0)
const EV_FORWARD  = Int8(1)
const EV_REVERSE  = Int8(2)

struct Event
    ii   :: Int32
    jj   :: Int32
    kind :: Int8
end

struct CellList
    cells     :: Vector{Vector{Int}}   # cell → particle indices
    cell_w    :: Float64
    cell_h    :: Float64
    ncols     :: Int
    nrows     :: Int
end
