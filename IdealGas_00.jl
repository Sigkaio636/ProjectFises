"""
Ideal Gas Simulation in Julia
==============================
- Toroidal (periodic) boundary conditions
- Particle struct with position and velocity
- Elastic hard-sphere collisions conserving momentum and kinetic energy
- Visualization via Plots.jl (GR backend):
    • snapshot frames saved as PNG
    • animation saved as GIF
    • thermodynamic time-series (KE, T, |p|) saved as PNG
    • final Maxwell-Boltzmann speed distribution saved as PNG

Dependencies (install once):
    julia> using Pkg
    julia> Pkg.add(["Plots", "Statistics", "Printf", "LinearAlgebra"])
"""

module IdealGas

using LinearAlgebra
using Statistics
using Printf
using Plots
gr()                          # GR backend — fast, no display needed, saves files

# ─────────────────────────────────────────────────
# Particle type
# ─────────────────────────────────────────────────
mutable struct Particle
    pos::Vector{Float64}   # [x, y]
    vel::Vector{Float64}   # [vx, vy]
    mass::Float64
    radius::Float64
end

Particle(pos, vel; mass=1.0, radius=0.5) = Particle(pos, vel, mass, radius)

kinetic_energy(p::Particle) = 0.5 * p.mass * dot(p.vel, p.vel)
speed(p::Particle)          = norm(p.vel)

# ─────────────────────────────────────────────────
# Simulation box — toroidal boundary conditions
# ─────────────────────────────────────────────────
struct Box
    width::Float64
    height::Float64
end

function wrap!(p::Particle, box::Box)
    p.pos[1] = mod(p.pos[1], box.width)
    p.pos[2] = mod(p.pos[2], box.height)
end

# Minimum-image displacement b − a on a torus
function min_image_delta(a::Particle, b::Particle, box::Box)
    dx = b.pos[1] - a.pos[1]
    dy = b.pos[2] - a.pos[2]
    dx -= box.width  * round(dx / box.width)
    dy -= box.height * round(dy / box.height)
    return dx, dy
end

# ─────────────────────────────────────────────────
# Elastic collision (positional correction + impulse)
# ─────────────────────────────────────────────────
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

# ─────────────────────────────────────────────────
# Integration step
# ─────────────────────────────────────────────────
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

# ─────────────────────────────────────────────────
# Initialization
# ─────────────────────────────────────────────────
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

# ─────────────────────────────────────────────────
# Diagnostics
# ─────────────────────────────────────────────────
function diagnostics(particles::Vector{Particle})
    KE  = sum(kinetic_energy(p) for p in particles)
    px  = sum(p.mass * p.vel[1] for p in particles)
    py  = sum(p.mass * p.vel[2] for p in particles)
    T   = KE / length(particles)
    return (; KE, px, py, T)
end

# ─────────────────────────────────────────────────
# Visualization helpers
# ─────────────────────────────────────────────────

"""
Speed → RGB colour on a blue–white–red ramp (cold = blue, fast = red).
"""
function speed_color(v::Float64, vmax::Float64)
    t = clamp(v / vmax, 0.0, 1.0)
    r = t
    g = 1.0 - abs(2t - 1.0)
    b = 1.0 - t
    return RGB(r, g, b)
end

"""
Render one snapshot of the particle positions.
Draws ghost copies at the opposite torus edges for a seamless look.
Returns a Plots.Plot — does NOT save or display by itself.
"""
function snapshot(particles::Vector{Particle}, box::Box, step_idx::Int)
    speeds = speed.(particles)
    vmax   = max(maximum(speeds), 1e-6)

    plt = plot(
        xlims        = (0, box.width),
        ylims        = (0, box.height),
        aspect_ratio = :equal,
        legend       = false,
        grid         = false,
        framestyle   = :box,
        background_color_inside = :black,
        background_color        = :black,
        foreground_color        = :white,
        title         = @sprintf("step %d", step_idx),
        titlefontcolor = :white,
        titlefontsize  = 10,
        tickfontcolor  = :white,
        size           = (520, 520),
    )

    for p in particles
        c = speed_color(speed(p), vmax)
        r = p.radius

        # Main position + torus ghost copies near edges
        draw_x = [p.pos[1]]
        draw_y = [p.pos[2]]
        p.pos[1] < r             && push!(draw_x, p.pos[1] + box.width)
        p.pos[1] > box.width - r && push!(draw_x, p.pos[1] - box.width)
        p.pos[2] < r             && push!(draw_y, p.pos[2] + box.height)
        p.pos[2] > box.height- r && push!(draw_y, p.pos[2] - box.height)

        for x in draw_x, y in draw_y
            scatter!(plt, [x], [y];
                markersize        = r * 0.52,
                markercolor       = c,
                markerstrokewidth = 0)
        end

        # Velocity arrow
        plot!(plt,
            [p.pos[1], p.pos[1] + p.vel[1] * 1.5],
            [p.pos[2], p.pos[2] + p.vel[2] * 1.5];
            lc=c, lw=0.8, alpha=0.55)
    end

    return plt
end

"""
Three-panel thermodynamics time-series plot.
Saves to `path` and returns the plot object.
"""
function plot_thermodynamics(times, KE_hist, T_hist, p_hist; path="thermodynamics.png")
    plt = plot(layout=(3,1), size=(700, 580),
        background_color = :white,
        left_margin  = 8Plots.mm,
        bottom_margin = 4Plots.mm)

    plot!(plt[1], times, KE_hist;
        lc=:steelblue, lw=1.5, legend=:topright,
        label="KE", ylabel="Total KE",
        title="Thermodynamic observables")
    hline!(plt[1], [mean(KE_hist)]; lc=:red, ls=:dash, lw=1, label="⟨KE⟩")

    plot!(plt[2], times, T_hist;
        lc=:darkorange, lw=1.5, legend=:topright,
        label="T", ylabel="Temperature")
    hline!(plt[2], [mean(T_hist)]; lc=:red, ls=:dash, lw=1, label="⟨T⟩")

    plot!(plt[3], times, p_hist;
        lc=:mediumseagreen, lw=1.2, legend=false,
        xlabel="Time", ylabel="|p| total",
        label="|p|")

    savefig(plt, path)
    println("  Saved thermo     → $path")
    return plt
end

"""
Speed distribution histogram overlaid with the theoretical
2-D Maxwell-Boltzmann distribution P(v) = (m/kT) v exp(−mv²/2kT).
Saves to `path` and returns the plot object.
"""
function plot_speed_distribution(particles::Vector{Particle};
                                  T_eq=1.0, mass=1.0, path="speed_dist.png")
    speeds  = speed.(particles)
    v_range = range(0, maximum(speeds)*1.4, length=300)
    mb      = (mass ./ T_eq) .* v_range .* exp.(-(mass .* v_range.^2) ./ (2T_eq))

    plt = plot(size=(600, 380),
        background_color = :white,
        left_margin  = 6Plots.mm,
        bottom_margin = 6Plots.mm)

    histogram!(plt, speeds;
        normalize  = :pdf,
        bins       = 18,
        fillcolor  = :steelblue,
        fillalpha  = 0.55,
        label      = "Simulation",
        xlabel     = "Speed  |v|",
        ylabel     = "P(v)",
        title      = "Speed distribution")

    plot!(plt, v_range, mb;
        lc=:red, lw=2.2,
        label=@sprintf("Maxwell-Boltzmann (T = %.3f)", T_eq))

    savefig(plt, path)
    println("  Saved speed dist → $path")
    return plt
end

# ─────────────────────────────────────────────────
# Main driver
# ─────────────────────────────────────────────────
"""
Run the full simulation and write all output files to `out_dir`:

  animation.gif       — particle movie (one frame every `save_every` steps)
  thermodynamics.png  — KE, temperature, and |p| vs time
  speed_dist.png      — final speed histogram vs Maxwell-Boltzmann
  snapshot_XXXX.png   — individual frames (only if save_frames=true)
"""
function run(;
    N           = 50,
    Lx          = 100.0,
    Ly          = 100.0,
    T_init      = 2.0,
    radius      = 2.0,
    mass        = 1.0,
    dt          = 0.2,
    n_steps     = 600,
    save_every  = 10,       # capture one animation frame every N steps
    out_dir     = "ideal_gas_output",
    save_frames = false,    # also write individual snapshot PNGs
    print_every = 100,
)
    isdir(out_dir) || mkpath(out_dir)

    box       = Box(Lx, Ly)
    particles = init_particles(N, box; T=T_init, mass=mass, radius=radius)

    println("="^60)
    println("  Ideal Gas Simulation — Julia + Plots.jl")
    println("="^60)
    @printf "  N=%d  box=%.0fx%.0f  r=%.2f  dt=%.3f  steps=%d\n\n" N Lx Ly radius dt n_steps

    # ── History buffers ───────────────────────────────────────────────────
    times   = Float64[]
    KE_hist = Float64[]
    T_hist  = Float64[]
    p_hist  = Float64[]

    # ── Animation ─────────────────────────────────────────────────────────
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

    # ── Save animation ────────────────────────────────────────────────────
    anim_path = joinpath(out_dir, "animation.gif")
    gif(anim, anim_path; fps=20)
    println("  Saved animation  → $anim_path")

    # ── Save thermodynamics plot ──────────────────────────────────────────
    plot_thermodynamics(times, KE_hist, T_hist, p_hist;
        path=joinpath(out_dir, "thermodynamics.png"))

    # ── Save speed distribution ───────────────────────────────────────────
    T_eq = mean(T_hist[max(1, end÷2):end])   # equilibrium temperature
    plot_speed_distribution(particles; T_eq=T_eq, mass=mass,
        path=joinpath(out_dir, "speed_dist.png"))

    println("\n  All outputs in: $(abspath(out_dir))/")
    println("="^60)
    return particles, box
end

end  # module IdealGas

# ─────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────
if abspath(PROGRAM_FILE) == @__FILE__
    IdealGas.run(
        N           = 100,
        n_steps     = 1200,
        save_every  = 5,
        radius      = 3.0,
        out_dir     = "ideal_gas_output",
        save_frames = false,
    )
end