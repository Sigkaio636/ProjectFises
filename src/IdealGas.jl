# ─────────────────────────────────────────────────
# src/IdealGas.jl
# ─────────────────────────────────────────────────

module IdealGas

using LinearAlgebra
using Statistics
using Printf
using Plots
gr()

# --- Inclusión de los submódulos de manera ordenada ---
include("Types.jl")
include("Physics.jl")
include("Visualization.jl")

function init_particles(N::Int, box::Box; T=1.0, mass=1.0, radius=0.5)
    particles = Particle[]
    cols = ceil(Int, sqrt(N))
    rows = ceil(Int, N / cols)
    dx   = box.width  / cols
    dy   = box.height / rows
    sv   = sqrt(T / mass)
    for k in 1:N
        col = (k-1) % cols
        row = (k-1) ÷ cols
        push!(particles, Particle(
            [(col + 0.5)*dx, (row + 0.5)*dy],
            [sv * randn(),   sv * randn()];
            mass=mass, radius=radius))
    end
    mean_vx = mean(p.vel[1] for p in particles)
    mean_vy = mean(p.vel[2] for p in particles)
    for p in particles
        p.vel[1] -= mean_vx
        p.vel[2] -= mean_vy
    end
    return particles
end

function run(;
    N           = 50,
    Lx          = 100.0,
    Ly          = 100.0,
    T_init      = 2.0,
    radius      = 2.0,
    mass        = 1.0,
    dt          = 0.2,
    n_steps     = 600,
    save_every  = 10,       
    out_dir     = "ideal_gas_output",
    save_frames = false,    
    print_every = 100,
)
    isdir(out_dir) || mkpath(out_dir)

    box       = Box(Lx, Ly)
    particles = init_particles(N, box; T=T_init, mass=mass, radius=radius)

    println("="^60)
    println("  Ideal Gas Simulation — Julia + Plots.jl")
    println("="^60)
    @printf "  N=%d  box=%.0fx%.0f  r=%.2f  dt=%.3f  steps=%d\n\n" N Lx Ly radius dt n_steps

    times   = Float64[]
    KE_hist = Float64[]
    T_hist  = Float64[]
    p_hist  = Float64[]

    anim        = Animation()
    frame_count = 0

    function record!(s)
        d = diagnostics(particles)
        push!(times,   s * dt)
        push!(KE_hist, d.KE)
        push!(T_hist,  d.T)
        push!(p_hist,  sqrt(d.px^2 + d.py^2))
        if s % save_every == 0
            frame_count += 1
            plt = snapshot(particles, box, s)
            frame(anim, plt)
            if save_frames
                p = joinpath(out_dir, @sprintf("snapshot_%04d.png", s))
                savefig(plt, p)
            end
        end
    end

    record!(0)
    d0 = diagnostics(particles)
    @printf "  Step %5d | KE=%8.3f | T=%.4f | |p|=%.2e\n" 0 d0.KE d0.T sqrt(d0.px^2 + d0.py^2)

    total_collisions = 0
    for s in 1:n_steps
        ncoll = step!(particles, box, dt)
        total_collisions += ncoll
        record!(s)
        if s % print_every == 0
            d = diagnostics(particles)
            @printf "  Step %5d | KE=%8.3f | T=%.4f | |p|=%.2e | coll=%d\n" s d.KE d.T sqrt(d.px^2 + d.py^2) ncoll
        end
    end

    println("\n  Total collisions : $total_collisions")
    println("  Frames captured  : $frame_count")

    anim_path = joinpath(out_dir, "animation.gif")
    gif(anim, anim_path; fps=20)
    println("  Saved animation  → $anim_path")

    plot_thermodynamics(times, KE_hist, T_hist, p_hist;
        path=joinpath(out_dir, "thermodynamics.png"))

    T_eq = mean(T_hist[max(1, end÷2):end])   
    plot_speed_distribution(particles; T_eq=T_eq, mass=mass,
        path=joinpath(out_dir, "speed_dist.png"))

    println("\n  All outputs in: $(abspath(out_dir))/")
    println("="^60)
    return particles, box
end

end  # module IdealGas