using Plots
include("plotConventions.jl")

#Plot the movement of the DSF Maxima
if true
    # Position of the maxima
    α_list=[0.0,0.2,0.45,0.7,1.0]
    pos_dyn=[[(1,0.13)],[(1,0.13)],[(0.5,1),(1.5,1)],[(0.65,0.73),(1.35,0.73)],[(0.65,0.67),(1.35,0.67)]]
    pos_vmc=[[(1,0)],[(1,0)],[(0.57,1),(1.43,1)],[(0.58,0.67),(1.42,0.67)],[(0.5,0.42),(1.5,0.42)]]
        
    # Assign a color per α
    colors = palette(:viridis, length(α_list))

    # Make the base plot (adjust to your heatmap or other data plot)
    plt = plot(title="DSF Maxima Positions", xlabel="k/π", ylabel="ω/J₁", legend=true, xlim=(0,2), ylim=(0,3.5))

    # Plot maxima as scatter points
    for (i, pos) in enumerate(pos_dyn)
        x_vals = [p[1] for p in pos]
        y_vals = [p[2] for p in pos]
        scatter!(plt, x_vals, y_vals, color=colors[i], label="α = $(α_list[i])", markersize=6)
    end
    for (i, pos) in enumerate(pos_vmc)
        x_vals = [p[1] for p in pos]
        y_vals = [p[2] for p in pos]
        scatter!(plt, x_vals, y_vals, color=colors[i], marker=:cross, markersize=6, label=false)
    end
    display(plt)
    #savefig("Images/DSF_pos_max.png")
end