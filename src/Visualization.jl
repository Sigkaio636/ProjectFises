# ─────────────────────────────────────────────────
# src/Visualization.jl
# ─────────────────────────────────────────────────

BURNOUT = 300

function snapshot(particles::Vector{Particle}, box::Box, step_idx::Int,
                  total_react_direct::Int, total_react_reverse::Int)
    
    n = length(particles)

    # Pre-allocate coordinate arrays for the 3 species
    xs = [Float64[] for _ in 1:3]
    ys = [Float64[] for _ in 1:3]
    
    # Pre-allocate NaN-separated arrays for velocity arrows
    vx_lines = [Float64[] for _ in 1:3]
    vy_lines = [Float64[] for _ in 1:3]

    # Use sizehint! to prevent dynamic array resizing
    for i in 1:3
        sizehint!(xs[i], n)
        sizehint!(ys[i], n)
        sizehint!(vx_lines[i], n * 3) # 3 elements per line (start, end, NaN)
        sizehint!(vy_lines[i], n * 3)
    end

    # Collect all positions, ghosts, and velocities in one pass
    @inbounds for p in particles
        sp = p.species
        r  = SPECIES[sp].radius
        px, py = p.x, p.y

        # ── OPTIMIZATION 2: Calculate ghosts without array allocations ──
        # Use NaN to represent "no ghost needed"
        gx = px < r ? px + box.width : (px > box.width - r ? px - box.width : NaN)
        gy = py < r ? py + box.height : (py > box.height - r ? py - box.height : NaN)

        # Always push the primary particle
        push!(xs[sp], px)
        push!(ys[sp], py)

        # Push edge ghosts (if they exist)
        if !isnan(gx)
            push!(xs[sp], gx)
            push!(ys[sp], py)
        end
        if !isnan(gy)
            push!(xs[sp], px)
            push!(ys[sp], gy)
        end
        # Push corner ghosts (if both X and Y wrap)
        if !isnan(gx) && !isnan(gy)
            push!(xs[sp], gx)
            push!(ys[sp], gy)
        end

        # ── Group velocity arrows ──
        push!(vx_lines[sp], px, px + p.vx * 1.2, NaN)
        push!(vy_lines[sp], py, py + p.vy * 1.2, NaN)
    end

    sp_labels = [
        "H₂O  (m=$(SPECIES[1].mass), r=$(SPECIES[1].radius), e=$(SPECIES[1].energy))",
        "H₃O⁺ (m=$(SPECIES[2].mass), r=$(SPECIES[2].radius), e=$(SPECIES[2].energy))",
        "OH⁻  (m=$(SPECIES[3].mass), r=$(SPECIES[3].radius), e=$(SPECIES[3].energy))",
    ]

    # Initialize the plot with all the styling
    plt = plot(
        xlims=(0, box.width), ylims=(0, box.height),
        aspect_ratio=:equal, legend=:topright,
        grid=false, framestyle=:box,
        background_color_inside=:white, background_color=:white, foreground_color=:black,
        title=@sprintf("step %d  |  rxn→: %d  |  rxn←: %d", step_idx, total_react_direct, total_react_reverse),
        titlefontsize=9, size=(580, 580)
    )

    # Plot everything in exactly 6 calls (3 for arrows, 3 for particles)
    for sp in 1:3
        # Draw velocity arrows first so they render underneath the particles
        if !isempty(vx_lines[sp])
            plot!(plt, vx_lines[sp], vy_lines[sp];
                  lc=SPECIES[sp].color, lw=0.7, alpha=0.45, label=false)
        end

        # Draw all particles and ghosts for this species
        if !isempty(xs[sp])
            scatter!(plt, xs[sp], ys[sp];
                markercolor=SPECIES[sp].color,
                markersize=SPECIES[sp].radius * 0.55,
                markerstrokewidth=0,
                label=sp_labels[sp])
        end
    end

    return plt
end

function plot_thermodynamics(times, KE_hist, T_hist, p_hist;
                              path="thermodynamics.png")
    plt = plot(layout=(3,1), size=(750, 680),
        background_color = :white,
        left_margin      = 8Plots.mm,
        bottom_margin    = 4Plots.mm)

    plot!(plt[1], times, KE_hist;
        lc=:steelblue, lw=1.5, label="KE", ylabel="Total KE",
        title="Thermodynamic observables — reactive ideal gas")
    mean_KE = mean(KE_hist[BURNOUT:end])
    std_KE = std(KE_hist[BURNOUT:end])
    hline!(plt[1], [mean_KE]; lc=:red, ls=:dash, lw=1, label=@sprintf("⟨KE⟩ = %.2f ± %.2f", mean_KE, std_KE))
    vline!(plt[1], [BURNOUT]; lc=:gray16, ls=:dot, lw=1, label="Burnout end")

    plot!(plt[2], times, T_hist;
        lc=:darkorange, lw=1.5, label="T", ylabel="Temperature")
    mean_T = mean(T_hist[BURNOUT:end])
    std_T = std(T_hist[BURNOUT:end])
    hline!(plt[2], [mean_T]; lc=:red, ls=:dash, lw=1, label=@sprintf("⟨T⟩ = %.2f ± %.2f", mean_T, std_T))
    vline!(plt[2], [BURNOUT]; lc=:gray16, ls=:dot, lw=1, label="Burnout end")

    plot!(plt[3], times, p_hist;
        lc=:mediumseagreen, lw=1.2, ylabel="|p| total", label="|p|")

    savefig(plt, path)
    println("  Saved thermo     → $path")
    return plt
end

function plot_speed_distribution(particles::Vector{Particle};
                                  T_eq=1.0, path="speed_dist.png")
    sp_colors = [SPECIES[k].color for k in 1:3]
    sp_names  = [SPECIES[k].name  for k in 1:3]

    plt = plot(size=(700, 420),
        background_color = :white,
        left_margin      = 6Plots.mm,
        bottom_margin    = 6Plots.mm,
        title            = "Speed distributions by species")

    for sp in 1:3
        group = filter(p -> p.species == sp, particles)
        isempty(group) && continue
        speeds  = speed.(group)                          # uses scalar fields via speed()
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

function plot_population_equilibrium(times,
        N_H2O_hist, N_H3O_hist, N_OH_hist,
        KE_H2O_hist, KE_H3O_hist, KE_OH_hist, Kc_hist;
        path="population.png")

    plt = plot(layout=(3,1), size=(750, 680),
        background_color = :white,
        left_margin      = 8Plots.mm,
        bottom_margin    = 4Plots.mm)

    mean_H2O = mean(N_H2O_hist[BURNOUT:end]); std_H2O = std(N_H2O_hist[BURNOUT:end])
    plot!(plt[1], times, N_H2O_hist; lc=:mediumseagreen, lw=1.5, label=@sprintf("N_H₂O = %.2f ± %.2f", mean_H2O, std_H2O),
        ylabel="Population", xlabel="Time", legend=:topright)
    mean_H3O = mean(N_H3O_hist[BURNOUT:end]); std_H3O = std(N_H3O_hist[BURNOUT:end])
    plot!(plt[1], times, N_H3O_hist; lc=:tomato,         lw=1.5, label=@sprintf("N_H₃O⁺ = %.2f ± %.2f", mean_H3O, std_H3O))
    mean_OH = mean(N_OH_hist[BURNOUT:end]); std_OH = std(N_OH_hist[BURNOUT:end])
    plot!(plt[1], times, N_OH_hist;  lc=:cornflowerblue, lw=1.5, label=@sprintf("N_OH⁻ = %.2f ± %.2f", mean_H3O, std_H3O))
    vline!(plt[1], [BURNOUT]; lc=:gray16, ls=:dot, lw=1, label="Burnout end")

    plot!(plt[2], times, KE_H2O_hist; lc=:mediumseagreen, lw=1.5, label="KE_H₂O",
        ylabel="KE per species", xlabel="Time", legend=:topright)
    plot!(plt[2], times, KE_H3O_hist; lc=:tomato,         lw=1.5, label="KE_H₃O⁺")
    plot!(plt[2], times, KE_OH_hist;  lc=:cornflowerblue, lw=1.5, label="KE_OH⁻")

    plot!(plt[3], times, Kc_hist; lc=:darkgoldenrod1, lw=1.5, label="Kc",
        ylabel="[H₃O⁺][OH⁻]/[H₂O]²", xlabel="Time", legend=:topright)

    savefig(plt, path)
    println("  Saved population → $path")
    return plt
end

function plot_kc_line(T_hist, Kc_hist;
                              path="kc_line.png")
    plt = plot(layout=(1,1), size=(750, 680),
        background_color = :white,
        left_margin      = 8Plots.mm,
        bottom_margin    = 4Plots.mm)

    lnKc  = log.(Kc_hist)
    inv_T = 1.0 ./ T_hist

    plot!(plt[1], inv_T[BURNOUT:end], lnKc[BURNOUT:end];
        lc=:steelblue, lw=1.5, xlabel="beta", ylabel="lnKc", marker='o', markersize=2,
        title="ln Kc = -Deps beta")

    savefig(plt, path)

    println("  Saved population → $path")
    return plt
end