# ─────────────────────────────────────────────────
# src/IdealGas.jl
# ─────────────────────────────────────────────────

module IdealGas

using LinearAlgebra
using Statistics
using Printf
using Plots
using Base.Threads: @threads, nthreads, threadid
gr()

# --- Inclusión de archivos internos ---
include("Types.jl")
include("Physics.jl")
include("Visualization.jl")

function init_particles(N_H2O::Int, N_H3O::Int, N_OH::Int, box::Box; T_init=2.0)
    N = N_H2O + N_H3O + N_OH
    species_list = vcat(fill(1, N_H2O), fill(2, N_H3O), fill(3, N_OH))

    cols = ceil(Int, sqrt(N))
    rows = ceil(Int, N / cols)
    dx   = box.width  / cols
    dy   = box.height / rows

    particles = Particle[]
    for k in 1:N
        sp  = species_list[k]
        m   = SPECIES[sp].mass
        sv  = sqrt(2*T_init / m)
        th_s = rand()*2*pi
        col = (k-1) % cols
        row = (k-1) ÷ cols
        push!(particles, Particle(
            (col + 0.5)*dx,          # x
            (row + 0.5)*dy,          # y
            sv * cos(th_s),          # vx
            sv * sin(th_s),          # vy
            Int32(sp)))
    end
    

    # Zero momentum — update field access:
    total_m  = sum(mass(p) for p in particles)
    mean_px  = sum(mass(p) * p.vx for p in particles) / total_m
    mean_py  = sum(mass(p) * p.vy for p in particles) / total_m
    KE_before = sum(0.5 * mass(p) * (p.vx^2 + p.vy^2) for p in particles)
    for p in particles;  p.vx -= mean_px;  p.vy -= mean_py;  end
    KE_after  = sum(0.5 * mass(p) * (p.vx^2 + p.vy^2) for p in particles)
    scale = sqrt(KE_before / KE_after)
    for p in particles;  p.vx *= scale;  p.vy *= scale;  end

    return particles
end

function run(;
    N_H2O         = 60,
    N_H3O         = 5,
    N_OH         = 5,
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
    particles = init_particles(N_H2O, N_H3O, N_OH, box; T_init=T_init)

    println("="^65)
    println("  Reactive Ideal Gas Simulation — Julia + Plots.jl")
    println("="^65)
    @printf "  N_H2O=%d  N_H3O=%d  N_OH=%d  p_react=%.3f\n" N_H2O N_H3O N_OH p_react
    @printf "  box=%.0fx%.0f  T_init=%.2f  dt=%.3f  steps=%d\n\n" Lx Ly T_init dt n_steps

    # History
    times   = Float64[]
    KE_hist = Float64[]
    T_hist  = Float64[]
    p_hist  = Float64[]
    N_H2O_hist = Int[]
    N_H3O_hist = Int[]
    N_OH_hist = Int[]
    KE_H2O_hist = Float64[]
    KE_H3O_hist = Float64[]
    KE_OH_hist = Float64[]
    Kc_hist = Float64[]

    anim            = Animation()
    frame_count     = 0
    total_elastic   = 0
    total_react_direct = 0
    total_react_reverse = 0

    function record!(s)
        d = diagnostics(particles)
        push!(times,   s * dt)
        push!(KE_hist, d.KE)
        push!(T_hist,  d.T)
        push!(p_hist,  sqrt(d.px^2 + d.py^2))
        push!(N_H2O_hist, d.N_H2O)
        push!(N_H3O_hist, d.N_H3O)
        push!(N_OH_hist, d.N_OH)
        push!(KE_H2O_hist, d.KE_H2O)
        push!(KE_H3O_hist, d.KE_H3O)
        push!(KE_OH_hist, d.KE_OH)
        push!(Kc_hist, d.Kc)

        if s % save_every == 0
            frame_count += 1
            plt = snapshot(particles, box, s, total_react_direct, total_react_reverse)
            frame(anim, plt)
            if save_frames
                p2 = joinpath(out_dir, @sprintf("snapshot_%04d.png", s))
                savefig(plt, p2)
            end
        end
    end

    record!(0)
    d0 = diagnostics(particles)
    @printf "  Step %5d | KE=%8.3f | T=%.4f | NH2O=%d NH3O=%d NOH=%d\n" 0 d0.KE d0.T d0.N_H2O d0.N_H3O d0.N_OH

    for s in 1:n_steps
        ne, nrd, nrr = step!(particles, box, dt, p_react)
        total_elastic   += ne
        total_react_direct += nrd
        total_react_reverse += nrr
        record!(s)
        if s % print_every == 0
            d = diagnostics(particles)
            @printf "  Step %5d | KE=%8.3f | T=%.4f | NH2O=%d NH3O=%d NOH=%d | col=%d | rxn_dir=%d  rxn_rev=%d\n" s d.KE d.T d.N_H2O d.N_H3O d.N_OH total_elastic total_react_direct total_react_reverse
        end
    end

    println("\n  Total elastic collisions : $total_elastic")
    println("  Total direct reactions (2 H2O → H3O + OH): $total_react_direct")
    println("  Total reverse reactions (H3O + OH → 2 H2O): $total_react_reverse")
    println("  Frames captured          : $frame_count")

    anim_path = joinpath(out_dir, "animation.gif")
    gif(anim, anim_path; fps=20)
    println("  Saved animation  → $anim_path")

    plot_thermodynamics(times, KE_hist, T_hist, p_hist;
        path=joinpath(out_dir, "thermodynamics.png"))

    T_eq = mean(T_hist[max(1, end÷2):end])
    plot_speed_distribution(particles; T_eq=T_eq,
        path=joinpath(out_dir, "speed_dist.png"))

    plot_population_equilibrium(times, N_H2O_hist, N_H3O_hist, N_OH_hist, KE_H2O_hist, KE_H3O_hist, KE_OH_hist, Kc_hist;
        path=joinpath(out_dir, "population.png"))

    println("\n  All outputs in: $(abspath(out_dir))/")
    println("="^65)
    return particles, box
end

end  # module IdealGas