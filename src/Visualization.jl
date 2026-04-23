# ─────────────────────────────────────────────────
# src/Visualization.jl
# ─────────────────────────────────────────────────

function snapshot(particles::Vector{Particle}, box::Box, step_idx::Int, total_react_direct::Int, total_react_reverse::Int)
    sp_colors = [SPECIES[1].color, SPECIES[2].color, SPECIES[3].color]
    sp_labels = ["H2O (m=$(SPECIES[1].mass), r=$(SPECIES[1].radius), e=$(SPECIES[1].energy))",
                 "H3O+ (m=$(SPECIES[2].mass), r=$(SPECIES[2].radius), e=$(SPECIES[2].energy))",
                 "OH- (m=$(SPECIES[3].mass), r=$(SPECIES[3].radius), e=$(SPECIES[3].energy))"]

    plt = plot(
        xlims        = (0, box.width),
        ylims        = (0, box.height),
        aspect_ratio = :equal,
        legend       = :topright,
        grid         = false,
        framestyle   = :box,
        background_color_inside = :white, # = :black,
        background_color        = :white, # = :black,
        foreground_color        = :black, # = :white,
        title          = @sprintf("step %d  |  rxn dir: %d  |  rxn rev: %d", step_idx, total_react_direct, total_react_reverse),
        titlefontcolor = :black, # = :white,
        titlefontsize  = 9,
        tickfontcolor  = :black, # = :white,
        legendfontcolor = :black, # = :white,
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

function plot_thermodynamics(times, KE_hist, T_hist, p_hist;
                              path="thermodynamics.png")
    plt = plot(layout=(3,1), size=(750, 680),
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

function plot_population_equilibrium(times, N_H2O_hist, N_H3O_hist, N_OH_hist, KE_H2O_hist, KE_H3O_hist, KE_OH_hist, Kc;
     path="population.png")
    plt = plot(layout=(3,1), size=(750, 680),
        background_color = :white,
        left_margin  = 8Plots.mm,
        bottom_margin = 4Plots.mm)

    plot!(plt[1], times, N_H2O_hist; lc=:mediumseagreen,          lw=1.5, label="N_H2O",
        ylabel="Population", xlabel="Time", legend=:topright)
    plot!(plt[1], times, N_H3O_hist; lc=:tomato,  lw=1.5, label="N_H3O")
    plot!(plt[1], times, N_OH_hist; lc=:cornflowerblue,  lw=1.5, label="N_OH")

    plot!(plt[2], times, KE_H2O_hist; lc=:mediumseagreen,          lw=1.5, label="KE_H2O",
        ylabel="KE per species", xlabel="Time", legend=:topright)
    plot!(plt[2], times, KE_H3O_hist; lc=:tomato,  lw=1.5, label="KE_H3O")
    plot!(plt[2], times, KE_OH_hist; lc=:cornflowerblue,  lw=1.5, label="KE_OH")

    plot!(plt[3], times, Kc; lc=:mediumseagreen,          lw=1.5, label="Kc",
        ylabel="[H3O+][OH-]/[H2O]^2", xlabel="Time", legend=:topright)

    savefig(plt, path)
    println("  Saved population     → $path")
    return plt
end