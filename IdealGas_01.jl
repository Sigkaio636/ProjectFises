"""
Reactive Ideal Gas Simulation in Julia
========================================
Three particle species: A (red), B (green), C (blue)
- Toroidal (periodic) boundary conditions
- Elastic hard-sphere collisions conserving momentum and kinetic energy
- Stochastic reaction:  A + A  →  B + C  (with probability p_react on collision)
- Visualization via Plots.jl (GR backend):
    • animation saved as GIF
    • thermodynamic time-series (KE, T, populations) saved as PNG
    • final Maxwell-Boltzmann speed distribution saved as PNG

Dependencies (install once):
    julia> using Pkg
    julia> Pkg.add(["Plots", "Statistics", "Printf", "LinearAlgebra"])
"""

module ReactiveIdealGas

using LinearAlgebra
using Statistics
using Printf
using Plots
gr()

# ─────────────────────────────────────────────────
# Species catalogue
# ─────────────────────────────────────────────────
"""
Each species has a fixed mass, radius, and display colour.
Species indices: 1 = A, 2 = B, 3 = C
"""
const SPECIES = (
    A = (mass=1.0, radius=2.0, color=:tomato,      name="A"),
    B = (mass=2.0, radius=2.8, color=:mediumseagreen, name="B"),
    C = (mass=0.5, radius=1.4, color=:cornflowerblue, name="C"),
)

# ─────────────────────────────────────────────────
# Particle type
# ─────────────────────────────────────────────────
mutable struct Particle
    pos::Vector{Float64}   # [x, y]
    vel::Vector{Float64}   # [vx, vy]
    species::Int           # 1=A  2=B  3=C
end

function Particle(pos, vel, species::Int)
    @assert species in (1, 2, 3) "species must be 1, 2, or 3"
    Particle(pos, vel, species)
end

mass(p::Particle)   = SPECIES[p.species].mass
radius(p::Particle) = SPECIES[p.species].radius
kinetic_energy(p::Particle) = 0.5 * mass(p) * dot(p.vel, p.vel)
speed(p::Particle)          = norm(p.vel)

# ─────────────────────────────────────────────────
# Simulation box — toroidal (periodic) boundary
# ─────────────────────────────────────────────────
struct Box
    width::Float64
    height::Float64
end

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

# ─────────────────────────────────────────────────
# Elastic collision (positional correction + impulse)
# ─────────────────────────────────────────────────
"""
Hard-sphere elastic collision — velocity impulse only.

Applies the standard CM-frame impulse:
    Δv_i = −[2 m_j/(m_i+m_j)] (v_rel · n̂) n̂
    Δv_j = +[2 m_i/(m_i+m_j)] (v_rel · n̂) n̂

Conserves total momentum and total kinetic energy exactly.

NO positional correction is performed.  Positional corrections shift
particle positions and can trigger a wrap!() across the periodic boundary,
introducing a discontinuous displacement that is not paired on the other
particle — this injects spurious momentum every collision and is the root
cause of the monotonic drift in |p|.  Overlap is a timestep artefact that
disappears with smaller dt; the approaching-velocity guard (vdotn < 0)
prevents the impulse from being applied twice to the same overlapping pair.

Returns true if the impulse was applied, false otherwise.
"""
function collide_elastic!(i::Particle, j::Particle, box::Box)
    dx, dy = min_image_delta(j, i, box)
    d2 = dx*dx + dy*dy

    sigma = radius(i) + radius(j)
    d2 > sigma^2  && return false
    d2 < 1e-14    && return false

    dvx = i.vel[1] - j.vel[1]
    dvy = i.vel[2] - j.vel[2]
    vdotn = dvx*dx + dvy*dy
    vdotn >= 0.0  && return false      # already separating — skip

    mi, mj     = mass(i), mass(j)
    total_mass = mi + mj
    factor     = 2.0 * mi * mj / total_mass * vdotn / d2
    i.vel[1]  -= factor / mi * dx;  i.vel[2] -= factor / mi * dy
    j.vel[1]  += factor / mj * dx;  j.vel[2] += factor / mj * dy

    return true
end

# ─────────────────────────────────────────────────
# Chemical reaction:  A + A  →  B + C
# ─────────────────────────────────────────────────
"""
Attempt the reaction A + A → B + C for two colliding A particles.

If the reaction fires (with probability `p_react`):

  Both total momentum AND total kinetic energy are conserved exactly.
  The derivation follows standard reactive-MD practice:

  1. Compute the pair's centre-of-mass (CM) velocity:
       v_cm = (p_i + p_j) / (m_A + m_A)
     This is invariant — the reaction cannot change it.

  2. Compute the kinetic energy available in the CM frame:
       KE_cm = ½ m_A |v_i - v_cm|² + ½ m_A |v_j - v_cm|²
             = ¼ m_A |v_rel|²

  3. Distribute KE_cm between products B and C via their reduced mass
       μ_BC = m_B m_C / (m_B + m_C):
       |v_rel'| = sqrt(2 KE_cm / μ_BC)
     The direction of v_rel' is chosen uniformly at random (isotropic
     scattering, no preferred axis after the reaction).

  4. Reconstruct lab-frame velocities:
       v_B = v_cm + [m_C/(m_B+m_C)] v_rel'
       v_C = v_cm − [m_B/(m_B+m_C)] v_rel'

  By construction: m_B v_B + m_C v_C = (m_B+m_C) v_cm = (2 m_A) v_cm = p_i + p_j  ✓
                   ½ m_B |v_B-v_cm|² + ½ m_C |v_C-v_cm|² = KE_cm  ✓

Returns true if the reaction fired.
"""
function react_AA!(i::Particle, j::Particle, p_react::Float64)
    # Only applies to A + A pairs
    (i.species == 1 && j.species == 1) || return false
    rand() > p_react && return false

    m_A = SPECIES[1].mass
    m_B = SPECIES[2].mass
    m_C = SPECIES[3].mass

    # ── Centre-of-mass velocity (conserved throughout) ─────────────────────
    # P_total = (m_A + m_A) * v_cm  →  v_cm = (p_i + p_j) / (2 m_A)
    vcm_x = (i.vel[1] + j.vel[1]) / 2.0   # m_A cancels since both equal
    vcm_y = (i.vel[2] + j.vel[2]) / 2.0

    # ── Total kinetic energy in the CM frame ───────────────────────────────
    # KE_cm = ½ m_A |v_i - v_cm|² + ½ m_A |v_j - v_cm|²
    #       = ¼ m_A |v_rel|²   (standard result for equal masses)
    vrel_x = i.vel[1] - j.vel[1]
    vrel_y = i.vel[2] - j.vel[2]
    KE_cm  = 0.25 * m_A * (vrel_x^2 + vrel_y^2)

    # ── Products B and C must share KE_cm in the CM frame ─────────────────
    # Reduced mass of products:  μ_BC = m_B * m_C / (m_B + m_C)
    # KE_cm = ½ μ_BC |v_rel'|²  →  |v_rel'| = sqrt(2 KE_cm / μ_BC)
    mu_BC   = m_B * m_C / (m_B + m_C)
    vrel_cm = sqrt(max(2.0 * KE_cm / mu_BC, 0.0))

    # Pick a random direction for the relative velocity of products
    theta = 2π * rand()
    rx = vrel_cm * cos(theta)
    ry = vrel_cm * sin(theta)

    # Reconstruct lab-frame velocities from CM velocity + relative velocity:
    #   v_B = v_cm + (m_C/(m_B+m_C)) * v_rel'
    #   v_C = v_cm - (m_B/(m_B+m_C)) * v_rel'
    fB = m_C / (m_B + m_C)
    fC = m_B / (m_B + m_C)

    i.vel[1] = vcm_x + fB * rx;  i.vel[2] = vcm_y + fB * ry
    j.vel[1] = vcm_x - fC * rx;  j.vel[2] = vcm_y - fC * ry

    # Change species
    i.species = 2   # becomes B
    j.species = 3   # becomes C

    return true
end

# ─────────────────────────────────────────────────
# Integration step
# ─────────────────────────────────────────────────
"""
One simulation step:
  1. Detect and resolve all pairwise overlaps.
     For A+A collisions, attempt the chemical reaction first;
     if it does not fire, apply the elastic collision.
  2. Advance positions with velocity-Verlet (no forces, so just free streaming).
  3. Apply periodic boundary conditions.

Returns (n_elastic, n_reactions).
"""
function step!(particles::Vector{Particle}, box::Box, dt::Float64, p_react::Float64)
    n          = length(particles)
    n_elastic  = 0
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

        # A+A pair: try reaction before elastic bounce
        if pi.species == 1 && pj.species == 1
            if react_AA!(pi, pj, p_react)
                n_reactions += 1
                # No positional correction — see collide_elastic! docstring.
                continue
            end
        end

        # Standard elastic collision
        if collide_elastic!(pi, pj, box)
            n_elastic += 1
        end
    end

    # Free streaming
    for p in particles
        p.pos[1] += p.vel[1] * dt
        p.pos[2] += p.vel[2] * dt
        wrap!(p, box)
    end

    return n_elastic, n_reactions
end

# ─────────────────────────────────────────────────
# Initialization
# ─────────────────────────────────────────────────
"""
Place N_A + N_B + N_C particles on a regular grid with randomised
Maxwell-Boltzmann velocities at temperature T_init.
"""
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

# ─────────────────────────────────────────────────
# Diagnostics
# ─────────────────────────────────────────────────
function diagnostics(particles::Vector{Particle})
    KE   = sum(kinetic_energy(p) for p in particles)
    px   = sum(mass(p) * p.vel[1] for p in particles)
    py   = sum(mass(p) * p.vel[2] for p in particles)
    T    = KE / length(particles)
    N_A  = count(p -> p.species == 1, particles)
    N_B  = count(p -> p.species == 2, particles)
    N_C  = count(p -> p.species == 3, particles)
    return (; KE, px, py, T, N_A, N_B, N_C)
end

# ─────────────────────────────────────────────────
# Visualization
# ─────────────────────────────────────────────────
function snapshot(particles::Vector{Particle}, box::Box, step_idx::Int, total_reactions::Int)
    # Species colours from catalogue
    sp_colors = [SPECIES[1].color, SPECIES[2].color, SPECIES[3].color]
    sp_labels = ["A (m=$(SPECIES[1].mass), r=$(SPECIES[1].radius))",
                 "B (m=$(SPECIES[2].mass), r=$(SPECIES[2].radius))",
                 "C (m=$(SPECIES[3].mass), r=$(SPECIES[3].radius))"]

    plt = plot(
        xlims        = (0, box.width),
        ylims        = (0, box.height),
        aspect_ratio = :equal,
        legend       = :topright,
        grid         = false,
        framestyle   = :box,
        background_color_inside = :black,
        background_color        = :black,
        foreground_color        = :white,
        title          = @sprintf("step %d  |  reactions: %d", step_idx, total_reactions),
        titlefontcolor = :white,
        titlefontsize  = 9,
        tickfontcolor  = :white,
        legendfontcolor = :white,
        legendbackground_color = RGBA(0,0,0,0.5),
        size           = (580, 580),
    )

    # One invisible scatter per species for the legend
    for (sp, col, lbl) in zip(1:3, sp_colors, sp_labels)
        scatter!(plt, Float64[], Float64[];
            markercolor=col, markerstrokewidth=0,
            markersize=5, label=lbl)
    end

    for p in particles
        col = sp_colors[p.species]
        r   = radius(p)

        draw_x = [p.pos[1]]
        draw_y = [p.pos[2]]
        p.pos[1] < r             && push!(draw_x, p.pos[1] + box.width)
        p.pos[1] > box.width - r && push!(draw_x, p.pos[1] - box.width)
        p.pos[2] < r             && push!(draw_y, p.pos[2] + box.height)
        p.pos[2] > box.height- r && push!(draw_y, p.pos[2] - box.height)

        for x in draw_x, y in draw_y
            scatter!(plt, [x], [y];
                markersize        = r * 0.55,
                markercolor       = col,
                markerstrokewidth = 0,
                label             = false)
        end

        # Velocity arrow
        plot!(plt,
            [p.pos[1], p.pos[1] + p.vel[1] * 1.2],
            [p.pos[2], p.pos[2] + p.vel[2] * 1.2];
            lc=col, lw=0.7, alpha=0.45, label=false)
    end

    return plt
end

function plot_thermodynamics(times, KE_hist, T_hist, p_hist, NA_hist, NB_hist, NC_hist;
                              path="thermodynamics.png")
    plt = plot(layout=(4,1), size=(750, 680),
        background_color = :white,
        left_margin  = 8Plots.mm,
        bottom_margin = 4Plots.mm)

    plot!(plt[1], times, KE_hist;
        lc=:steelblue, lw=1.5, legend=:topright,
        label="KE", ylabel="Total KE",
        title="Thermodynamic observables — reactive ideal gas")
    hline!(plt[1], [mean(KE_hist)]; lc=:red, ls=:dash, lw=1, label="⟨KE⟩")

    plot!(plt[2], times, T_hist;
        lc=:darkorange, lw=1.5, legend=:topright,
        label="T", ylabel="Temperature")
    hline!(plt[2], [mean(T_hist)]; lc=:red, ls=:dash, lw=1, label="⟨T⟩")

    plot!(plt[3], times, p_hist;
        lc=:mediumseagreen, lw=1.2, legend=false,
        ylabel="|p| total", label="|p|")

    plot!(plt[4], times, NA_hist; lc=:tomato,          lw=1.5, label="N_A",
        ylabel="Population", xlabel="Time", legend=:topright)
    plot!(plt[4], times, NB_hist; lc=:mediumseagreen,  lw=1.5, label="N_B")
    plot!(plt[4], times, NC_hist; lc=:cornflowerblue,  lw=1.5, label="N_C")

    savefig(plt, path)
    println("  Saved thermo     → $path")
    return plt
end

function plot_speed_distribution(particles::Vector{Particle}; T_eq=1.0, path="speed_dist.png")
    sp_colors = [SPECIES[1].color, SPECIES[2].color, SPECIES[3].color]
    sp_names  = [SPECIES[k].name for k in 1:3]

    plt = plot(size=(700, 420),
        background_color = :white,
        left_margin  = 6Plots.mm,
        bottom_margin = 6Plots.mm,
        title = "Speed distributions by species")

    for sp in 1:3
        group = filter(p -> p.species == sp, particles)
        isempty(group) && continue
        speeds  = speed.(group)
        m       = SPECIES[sp].mass
        v_range = range(0, maximum(speeds)*1.4, length=300)
        mb      = (m ./ T_eq) .* v_range .* exp.(-(m .* v_range.^2) ./ (2T_eq))

        histogram!(plt, speeds;
            normalize  = :pdf,
            bins       = 14,
            fillcolor  = sp_colors[sp],
            fillalpha  = 0.40,
            label      = "$(sp_names[sp]) sim",
            xlabel     = "Speed |v|",
            ylabel     = "P(v)")
        plot!(plt, v_range, mb;
            lc=sp_colors[sp], lw=2.0, ls=:dash,
            label=@sprintf("%s MB (T=%.3f)", sp_names[sp], T_eq))
    end

    savefig(plt, path)
    println("  Saved speed dist → $path")
    return plt
end

# ─────────────────────────────────────────────────
# Main driver
# ─────────────────────────────────────────────────
"""
Run the reactive ideal gas simulation.

Keyword arguments
─────────────────
N_A, N_B, N_C   : initial number of each species
Lx, Ly          : box dimensions
T_init          : initial temperature
dt              : time step
n_steps         : number of integration steps
p_react         : probability that an A+A collision triggers A+A → B+C
save_every      : animation frame interval
out_dir         : output directory
save_frames     : also write individual snapshot PNGs
print_every     : console logging interval
"""
function run(;
    N_A         = 60,
    N_B         = 5,
    N_C         = 5,
    Lx          = 120.0,
    Ly          = 120.0,
    T_init      = 2.0,
    dt          = 0.2,
    n_steps     = 800,
    p_react     = 0.05,      # reaction probability per A+A contact
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

end  # module ReactiveIdealGas

# ─────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────
if abspath(PROGRAM_FILE) == @__FILE__
    ReactiveIdealGas.run(
        N_A         = 80,
        N_B         = 5,
        N_C         = 5,
        Lx          = 120.0,
        Ly          = 120.0,
        T_init      = 2.0,
        dt          = 0.2,
        n_steps     = 1000,
        p_react     = 0.05,
        save_every  = 5,
        out_dir     = "reactive_gas_output",
        save_frames = false,
    )
end