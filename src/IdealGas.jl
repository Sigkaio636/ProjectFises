# ─────────────────────────────────────────────────
# src/IdealGas.jl
# ─────────────────────────────────────────────────

module IdealGas

using LinearAlgebra
using Statistics
using Printf
using Plots
gr()

# --- Inclusión de archivos internos ---
include("Types.jl")
include("Physics.jl")
include("Visualization.jl")

function init_particles(N_A::Int, N_B::Int, N_C::Int, box::Box; T_init=2.0)
    N = N_A + N_B + N_C
    species_list = vcat(fill(1, N_A), fill(2, N_B), fill(3, N_C))

    cols = ceil(Int, sqrt(N))
    rows = ceil(Int, N / cols)
    dx   = box.width  / cols
    dy   = box.height / rows

    particles = Particle[]
    for k in 1:N
        sp  = species_list[k]
        m   = SPECIES[sp].mass
        sv  = sqrt(T_init / m)
        col = (k-1) % cols
        row = (k-1) ÷ cols
        push!(particles, Particle(
            [(col + 0.5)*dx, (row + 0.5)*dy],
            [sv * randn(),   sv * randn()],
            sp))
    end

    # Zero total momentum
    total_m  = sum(mass(p) for p in particles)
    mean_px  = sum(mass(p) * p.vel[1] for p in particles) / total_m
    mean_py  = sum(mass(p) * p.vel[2] for p in particles) / total_m
    for p in particles
        p.vel[1] -= mean_px
        p.vel[2] -= mean_py
    end
    return particles
end

function run(;
    N_A         = 60,
    N_B         = 5,
    N_C         = 5,
    Lx          = 120.0,
    Ly          = 120.0,
    T_init      = 2.0,
    dt          = 0.2,
    n_steps     = 800,
    p_react     = 0.05,
    save_every  = 10,
    out_dir     = "reactive_gas_output",
    save_frames = false,
    print_every = 100,
)
    isdir(out_dir) || mkpath(out_dir)

    box       = Box(Lx, Ly)
    particles = init_particles(N_A, N_B, N_C, box; T_init=T_init)

    println("="^65)
    println("  Reactive Ideal Gas Simulation — Julia + Plots.jl")
    println("="^65)
    @printf "  N_A=%d  N_B=%d  N_C=%d  p_react=%.3f\n" N_A N_B N_C p_react
    @printf "  box=%.0fx%.0f  T_init=%.2f  dt=%.3f  steps=%d\n\n" Lx Ly T_init dt n_steps

    # History
    times   = Float64[]
    KE_hist = Float64[]
    T_hist  = Float64[]
    p_hist  = Float64[]
    NA_hist = Int[]
    NB_hist = Int[]
    NC_hist = Int[]

    anim            = Animation()
    frame_count     = 0
    total_elastic   = 0
    total_reactions = 0

    function record!(s)
        d = diagnostics(particles)
        push!(times,   s * dt)
        push!(KE_hist, d.KE)
        push!(T_hist,  d.T)
        push!(p_hist,  sqrt(d.px^2 + d.py^2))
        push!(NA_hist, d.N_A)
        push!(NB_hist, d.N_B)
        push!(NC_hist, d.N_C)
        if s % save_every == 0
            frame_count += 1
            plt = snapshot(particles, box, s, total_reactions)
            frame(anim, plt)
            if save_frames
                p2 = joinpath(out_dir, @sprintf("snapshot_%04d.png", s))
                savefig(plt, p2)
            end
        end
    end

    record!(0)
    d0 = diagnostics(particles)
    @printf "  Step %5d | KE=%8.3f | T=%.4f | NA=%d NB=%d NC=%d\n" 0 d0.KE d0.T d0.N_A d0.N_B d0.N_C

    for s in 1:n_steps
        ne, nr = step!(particles, box, dt, p_react)
        total_elastic   += ne
        total_reactions += nr
        record!(s)
        if s % print_every == 0
            d = diagnostics(particles)
            @printf "  Step %5d | KE=%8.3f | T=%.4f | NA=%d NB=%d NC=%d | rxn=%d\n" s d.KE d.T d.N_A d.N_B d.N_C total_reactions
        end
    end

    println("\n  Total elastic collisions : $total_elastic")
    println("  Total reactions (A+A→B+C): $total_reactions")
    println("  Frames captured          : $frame_count")

    anim_path = joinpath(out_dir, "animation.gif")
    gif(anim, anim_path; fps=20)
    println("  Saved animation  → $anim_path")

    plot_thermodynamics(times, KE_hist, T_hist, p_hist, NA_hist, NB_hist, NC_hist;
        path=joinpath(out_dir, "thermodynamics.png"))

    T_eq = mean(T_hist[max(1, end÷2):end])
    plot_speed_distribution(particles; T_eq=T_eq,
        path=joinpath(out_dir, "speed_dist.png"))

    println("\n  All outputs in: $(abspath(out_dir))/")
    println("="^65)
    return particles, box
end

end  # module IdealGas