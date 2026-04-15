# ─────────────────────────────────────────────────
# src/Visualization.jl
# ─────────────────────────────────────────────────

function speed_color(v::Float64, vmax::Float64)
    t = clamp(v / vmax, 0.0, 1.0)
    r = t
    g = 1.0 - abs(2t - 1.0)
    b = 1.0 - t
    return RGB(r, g, b)
end

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