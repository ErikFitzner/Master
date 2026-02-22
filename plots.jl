using Plots

include("plotConventions.jl")

#=
plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S(w=0)/J", size=(0.8*1.5*aps_width,0.8*aps_width),legend=:bottomright, xlims=(1.9,6.1), ylims=(0.25,1.05))
r_max = [1,2,3,4,5]
ys = [0.6786602297873916, 0.7039581360308462, 0.7166268745594179, 0.7151146533342323, 0.7220749266307802]
yl = [0.9594092957640764, 0.995147257137703, 0.9813112720447729, 0.9239983658010862, 0.940292695140516]
for r in r_max
    if r == 1
        Plots.scatter!(plt5, [r+1], [yl[r]], color=color_vec[r], markersize=6, label="ladder")
        Plots.scatter!(plt5, [r+1], [ys[r]], color=color_vec[r], markersize=7, marker=:cross, label="square lattice")
    else
        Plots.scatter!(plt5, [r+1], [yl[r]], color=color_vec[r], markersize=6, label=false)
        Plots.scatter!(plt5, [r+1], [ys[r]], color=color_vec[r], markersize=7, marker=:cross, label=false)
    end
end
# Plot shaded band using ribbon
xvals = [1.9, 6.1]
ymid = (0.98 + 0.96) / 2
yribbon = (0.98 - 0.96) / 2
Plots.plot!(plt5, xvals, fill(ymid, 2), ribbon=fill(yribbon, 2), fillalpha=0.2, color=color_vec[7], linestyle=:dash, label="[Rakovszky...2022]")
Plots.hline!([0.95], color=color_vec[8], linestyle=:dash, label="[Kloss...2018]", linewidth=0.8)
Plots.hline!([0.72], color=:black, linestyle=:dash, label=false, linewidth=0.8)
display(plt5)
# savefig(plt5,"Images/SpinDiffusion_XX_n.png")
=#

#=
plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J(x=0)",title="Spin Diffusion constant",size=(1.5*aps_width,aps_width),legend=:bottomright, xlims=(1.9,6.1))  # ,title="Spin Diffusion, $(lattice_type)", ylabel=L"D_S/J\;(T\to∞)"

Plots.scatter!(plt5, [2], [sqrt(pi*bc1/(4*2-2-1/(4*bc1)))], label="[Kopietz1993]", markersize=10, color=:orange)
Plots.scatter!(plt5, [2.0,3.0,4.0], [0.793,0.842,0.839]./2, label="[Morita1972]", markershape=:star, markersize=14, color=:green)

Plots.scatter!(plt5, [2], [sqrt(pi*bc1/(4*3-2-1/(4*bc1)))], label=false, markersize=10, color=:orange)
Plots.scatter!(plt5, [2.0,3.0,4.0], [0.591,0.620,0.618]./2, label=false, markershape=:star, markersize=14, color=:green)

Plots.scatter!(plt5, [2.0,3.0,4.0], [0.492,0.512,0.507]./2, label=false, markershape=:star, markersize=14, color=:green)

Plots.scatter!(plt5, [2,3,4,5,6], [0.39633272976014317, 0.44311346272378266, 0.41970335222826516, 0.4192399728370765, 0.41690808857603734], label="Dyn-HTE", alpha=0.7, color=:black, markersize=6)
Plots.scatter!(plt5, [2,3,4,5], [0.2954089751508137, 0.3228543656794131, 0.30676251447626796, 0.31454842340367045], label=false, alpha=0.7, color=:black, markersize=6)
Plots.scatter!(plt5, [2,3,4,5], [0.24579512472432255, 0.2652084044513992, 0.24762055795367333, 0.25681786636263026], label=false, alpha=0.7, color=:black, markersize=6)

Plots.hline!([0.42], color=:red, linestyle=:dash, label="square", linewidth=0.8)
Plots.hline!([0.619/2], color=:magenta, linestyle=:dash, label="simple cubic", linewidth=0.8)
Plots.hline!([0.509/2], color=:blue, linestyle=:dash, label="bcc", linewidth=0.8)
display(plt5)
#savefig(plt5,"Images/SpinDiffusion_n.png")
=#


#=
plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J(w=0)",size=(0.8*1.5*aps_width,0.8*aps_width),legend=:bottomright, xlims=(1.9,6.1), ylims=(0.25,1.05)) # title="Spin Diffusion"
Plots.scatter!(plt5, [2,3,4,5,6], [0.9594092957640764, 0.995147257137703, 0.9813112720447729, 0.9239983658010862, 0.940292695140516], label="ladder", alpha=0.7, color=color_vec, markersize=6)
Plots.scatter!(plt5, [2,3,4,5,6], [0.6786602297873916, 0.7039581360308462, 0.7166268745594179, 0.7151146533342323, 0.7220749266307802], label="square lattice", alpha=0.7, color=color_vec, markershape=:cross, markersize=7)
# Plot shaded band using ribbon
xvals = [1.9, 6.1]
ymid = (0.98 + 0.96) / 2
yribbon = (0.98 - 0.96) / 2
Plots.plot!(plt5, xvals, fill(ymid, 2), ribbon=fill(yribbon, 2), fillalpha=0.2, color=color_vec[7], linestyle=:dash, label="[Pollmann2022]")
Plots.hline!([0.95], color=color_vec[8], linestyle=:dash, label="[Reichman2018]", linewidth=0.8)
Plots.hline!([0.72], color=:black, linestyle=:dash, label=false, linewidth=0.8)
display(plt5)
#savefig(plt5,"Images/SpinDiffusion_XX_n.png")
=#
