using Plots
using DelimitedFiles
using Polynomials#, RationalFunctions
using RobustPade
using Statistics

include("../../plotConventions.jl")
include("../../LatticeGraphs.jl")
include("../../Embedding.jl")
include("../../ConvenienceFunctions.jl")

# predifine variables and load moments
x_vec_sh = collect(0.1:0.01:10)
x0_vec_sh = 1 ./ [3.0,1.8,1.2,0.95,0.8,0.7,0.6,0.5,0.43,0.38]
n_max = 12
r_max = 6

m_vec_lines = readlines("m_vec.txt")
m_vec_loaded = [parse.(Float64, split(line)) for line in m_vec_lines]

#m0_vec_lines = readlines("m0_vec.txt")
#m0_vec_loaded = [parse.(Float64, split(line)) for line in m0_vec_lines]

# turn coefficients into polynomial, save them and plot moments
if true
    plt_sh = plot(xlabel="x=J/T", ylabel="m_k")
    data_vecs = Vector{Polynomial}(undef, r_max + 1)
    for r in 0:r_max
        plot!(plt_sh,x_vec_sh,Polynomial(m_vec_loaded[r+1]).(x_vec_sh))
        data_vecs[r+1] = Polynomial(m_vec_loaded[r+1])
    end
    #display(plt_sh)
end

# introduce function for QFI and plot the series x- and u-Padés
if true
    x_poly = Polynomial([0, 1]) 
    f_Q_poly = 
    x_poly^2    * data_vecs[2] -
    x_poly^4/12 * data_vecs[3] +
    x_poly^6/120 * data_vecs[4] -
    (34 * x_poly^8 / 40320) * data_vecs[5] +
    (62 * x_poly^10 / 725760) * data_vecs[6] -
    (1382 * x_poly^12 / 159667200) * data_vecs[7]

    # plot
    plt_sh3 = plot(xlabel="J/T", ylabel="f_Q[k,T]",legend=:topleft, ylims=(0,1))
    plot!(plt_sh3,x_vec_sh,f_Q_poly.(x_vec_sh),color="gray",label="x-Series")

    # x-Pade
    plot!(plt_sh3,x_vec_sh,robustpade(f_Q_poly,5,6).(x_vec_sh),color="blue",alpha=0.7,label="x-Padé[5,6]")
    plot!(plt_sh3,x_vec_sh,robustpade(f_Q_poly,4,5).(x_vec_sh),color="red",alpha=0.7,label="x-Padé[4,5]")

    # u-Pade
    f=0.23  #f=0.23 for [4,5] and [5,6], 0.3 für [5,5],[6,6] (wähle f so, dass möglichst viele Padés übereinander liegen)
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
    u_vec = tanh.(f .* x_vec_sh)
    coeffs_x = coeffs(f_Q_poly)
    p_u = Polynomial(ufromx_mat*coeffs_x)


    # check for underflow/overflow (for small f)
    if false
        println("Min abs u_vec: ", minimum(abs.(u_vec)))
        println("Min abs ufromx_mat: ", minimum(abs.(ufromx_mat)))

        num_zeros_u = count(iszero, u_vec)
        println("Zeros in u_vec: ", num_zeros_u)
        num_zeros_mat = count(iszero, ufromx_mat)
        println("Zeros in ufromx_mat: ", num_zeros_mat)

        println(coeffs_x)
        println(ufromx_mat*coeffs_x)
    end

    #plot!(plt_sh3,x_vec_sh,robustpade(p_u,3,4).(u_vec),color="green",linestyle=:dash, label="u-Padé[3,4]")
    plot!(plt_sh3,x_vec_sh,robustpade(p_u,4,5).(u_vec),color="blue",linestyle=:dash,alpha=0.7,label="u-Padé[4,5] (f=$f)")
    plot!(plt_sh3,x_vec_sh,robustpade(p_u,5,6).(u_vec),color="red",linestyle=:dashdot,alpha=0.7,label="u-Padé[5,6] (f=$f)")
    #display(plt_sh3)
    #savefig(plt_sh3,"CaseStudy/Triangular_Lattice/Shannon_f.png")
end

# plot against 1/x and compare to results from [Shimokawa...2025]
if true
    inv_x_vec = 1.0 ./ x_vec_sh
    
    x_sh = collect(0.1:0.1:2.0)
    y_sh = [3.706,2.397,1.667,1.314,1.079,0.929,0.778,0.683,0.595,0.536,0.485,0.440,0.397,0.357,0.324,0.298,0.273,0.252,0.237,0.221] # N=30
    y_sh2 = [3.544,2.470,1.746,1.346,1.089,0.922,0.784,0.683,0.595,0.536,0.485,0.440,0.397,0.357,0.324,0.298,0.273,0.252,0.237,0.221] # N=24
    y_sh3 = [3.192,2.451,1.759,1.384,1.149,0.966,0.829,0.717,0.630,0.556,0.486,0.440,0.397,0.357,0.324,0.298,0.273,0.252,0.237,0.221] # N=18
    
    # Plot
    #plt_sh4 = plot(xlabel="T/J", ylabel="f_Q[k,T]",xlims=(0,2.5), ylims=(0,4))
    
    # log-log
    plt_sh4 = plot(xlabel="T/J", ylabel="f_Q[k,T]",xscale=:log10,yscale=:log10,xlims=(0.08,1), ylims=(0.5,10))

    #plot!(plt_sh4, inv_x_vec, f_Q_poly.(x_vec_sh),color="gray",label="x-Series")
    
    #plot!(plt_sh4,inv_x_vec,robustpade(f_Q_poly,5,6).(x_vec_sh),color="blue",label="x-Padé[5,6]")
    #plot!(plt_sh4,inv_x_vec,robustpade(f_Q_poly,4,5).(x_vec_sh),color="red",label="x-Padé[4,5]")
    
    plot!(plt_sh4,inv_x_vec,robustpade(p_u,5,6).(u_vec),color="blue",linestyle=:dash,alpha=0.7,label="u-Padé[5,6] (f=$f)")
    plot!(plt_sh4,inv_x_vec,robustpade(p_u,4,5).(u_vec),color="red",linestyle=:dashdot,alpha=0.7,label="u-Padé[4,5] (f=$f)")

    # Powerlaw
    # Model: y = A * x^α
    powerlaw(x, p) = p[1] .* x .^ p[2]  # p[1] = A, p[2] = α

    xdata = collect(1:0.01:2.5)
    ydata = robustpade(p_u,4,5).(tanh.(f .* xdata))

    p0 = [1.0, -1.0]  # initial guess [A, α]
    fit = curve_fit(powerlaw, xdata, ydata, p0)

    A, α = fit.param
    println("Fitted A = $A")
    println("Fitted alpha = $α")
    # Generate smooth x values for fit line
    xfine = range(0.4, stop=10, length=200)
    yfit = powerlaw(xfine, fit.param)

    # Plot fitted curve
    plot!(xfine.^(-1), yfit, label="Fit: y = $(round(A, digits=2))·x^($(round(α, digits=2)))", lw=2)

    #plot!(plt_sh4,x_sh,0.52*x_sh.^(-1),color=:gray,label="y~1/x")
    #pl(x) = 0.52*x^(1.06) #0.52+0.61*(x-1)
    #plot!(plt_sh4,inv_x_vec,pl.(x_vec_sh),color=:gray,label="y~1/x")

    plot!(plt_sh4,[NaN],[NaN],color=:white,label=" ",linealpha=0,marker=:none)
    plot!(plt_sh4,[NaN],[NaN],color=:white,label="[Shimokawa...2025]",linealpha=0,marker=:none)
    
    plot!(plt_sh4,x_sh,y_sh,color="orange",marker=:circle,markersize=6,alpha=0.3,label="N=30")
    plot!(plt_sh4,x_sh,y_sh2,color="orange",marker=:rect,markersize=6,alpha=0.3,label="N=24")
    plot!(plt_sh4,x_sh,y_sh3,color="orange",marker=:utriangle,markersize=6,alpha=0.3,label="N=18")
    display(plt_sh4)
    # savefig(plt_sh4,"CaseStudy/Triangular_Lattice/Shannon_x.png")

    #println(robustpade(p_u,4,5))
end

# ratio method to extract power law behavior 
if false
    n = range(5,12)
    r_n = [coeffs_x[i+1]/coeffs_x[i] for i in n]
    println(coeffs_x)

    n = 1 ./ n
    a = cor(n, r_n) * (std(r_n) / std(n))
    b = mean(r_n) - a * mean(n)
    println("Best fit line: y = $a * x + $b")
    
    x_c = 1/b
    theta = x_c*a+1
    println("x_c = $x_c, theta = $theta")

    # Plot
    plt_sh5 = plot(xlabel="1/n", ylabel="r_n")
    scatter!(plt_sh5,n,r_n)
    plot!(plt_sh5,n, a .* n .+ b, label="y=$(round(a, digits=2))*x + $(round(b, digits=2))", color=:red)
    display(plt_sh5)
end

# Dlog Padé
if false
    inv_x_vec = 1.0 ./ x_vec_sh
    f_der = Polynomials.derivative(f_Q_poly)
    log_der(x) = f_der(x) / f_Q_poly(x)
    plt_sh6 = plot(xlabel="T/J", ylabel="(f_Q'(x))/f_Q(x)",xlims=(0,2.5),ylims=(0,15))
    plot!(plt_sh6,inv_x_vec, log_der.(x_vec_sh), label="DLog")

    
    # x-Padé
    # Choose expansion point x_0 away from roots of denominator
    #println(roots(f_Q_poly))
    x_0 = 1.5

    plot!(plt_sh6,inv_x_vec,robustpade(log_der,5,6,x_0).(x_vec_sh),color="blue",alpha=0.7,label="x-Padé[5,6]")
    plot!(plt_sh6,inv_x_vec,robustpade(log_der,4,5,x_0).(x_vec_sh),color="red",alpha=0.7,label="x-Padé[4,5]")
    #display(plt_sh6)
end