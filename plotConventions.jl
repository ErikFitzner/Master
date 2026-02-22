using Plots, LaTeXStrings, ColorSchemes, Colors
pyplot()    # using gr plotting backend
Plots.default(titlefont=(10,), legendfontsize=8, background_color_legend=nothing, foreground_color_legend=nothing,
grid=false, guidefont=(9,), tickfont=(9,), framestyle=:box, linewidth=1,markerstrokewidth=0,markersize=3.5,markeralpha=0.9, dpi=200)
aps_width = 325 #in pixels, seems to work with pyplot backend  (243 # width of single-colum plot in points, it assumes internally dpi=100 (so we set this above), but this does not limit the resolution since we later save as pdf
color_vec = [palette(:tab10)[c%10+1] for c in 0:50] #
#color_vec = ["blue","green","red","cyan","brown","purple","magenta","orange","teal","aquamarine2","steelblue1","darkorchid","grey56","olive","blue","green","red","cyan","brown","purple","magenta","orange","teal","aquamarine2","steelblue1","darkorchid","grey56","olive"]
grey_vec = ["black","grey","darkblue","olive","brown","steelblue1"]
thermalCol4_vec = reverse(ColorSchemes.thermal[1:70:end])
thermalCol13_vec = reverse(ColorSchemes.thermal[1:20:end])
marker_vec = [:dot,:cross,:diamond,:dtriangle,:square,:star4,:utriangle]
linestyle_vec = [:solid, :dash, :dot, :dashdot, :solid, :dash, :dot, :dashdot]

plt_empty = Plots.plot([0,0],[0,0],label="")

function saturation_palette(N; h=240, v=1.0, s_min=0.2, s_max=1.0)
    ss = range(s_min, s_max; length=N)
    return [RGB(HSV(h, s, v)) for s in ss]
end

function addABC(plt,label);
    """ adds "label" on the top left of plot axes """
    xl,yl = xlims(plt),ylims(plt)
    annotate!(plt, xl[1]-0.16*(xl[2]-xl[1]), yl[2]-0.02*(yl[2]-yl[1]), (label,8,:left) )
end

# Exact functions for J1-J2 tetramer
function exact11(x, a) 
    result = (
    12 * (-2 + a) * (-1 + a) * (1 + a) * exp(x + 2a * x) +
    a * (
        4 * (-1 + a) * (1 + a) * (-5 + 2a) * exp(3x) +
        3 * (-2 + a) * (1 + a) * exp(2x) * (-3 + a * (-1 + x) - x) +
        5 * (-2 + a) * (-1 + a) * (-5 - a + 3 * (1 + a) * x)
    ) +
    12 * exp(x + a * x) * (
        -2 + (-1 + a) * a * (
            2 * (-1 + a) + (-2 + a) * (1 + a) * x
        )
    )
    ) / (
    24 * (-2 + a) * (-1 + a) * a * (1 + a) * (
        5 + exp(x) * (
            6 * exp(a * x) +
            exp(2a * x) +
            exp(x) * (3 + exp(x))
        )
    ) * x
    )
    return result
end
function exact12(x, a)
    numerator = 5 + 15*x + exp(2*x)*(3 - 8*exp(x) + 3*x)
    denominator = 24 * (5 + 3*exp(2*x) + exp(3*x) + 6*exp((1 + a)*x) + exp(x + 2*a*x)) * x
    return numerator / denominator
end
