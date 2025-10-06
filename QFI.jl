using JLD2, DelimitedFiles, SimpleWeightedGraphs, Plots

include("plotConventions.jl")
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 

### load graph evaluations
spin_length = 1/2
bc1 = spin_length*(spin_length+1)/3
n_max = 12

### prepare lattice
lattice_type = "Shastry-Sutherland" # "Shastry-Sutherland"

# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = true
j3 = false
j4 = false

L = 12

hte_lattice = getLattice(L,lattice_type,j1,j2,j3,j4);

### compute all correlations in the lattice (or load them)
function load_ciipDyn_mat_subst(a::Float64;b::Float64=0.0,c::Float64=0.0)
    fileName_c = "CaseStudy/$(lattice_type)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_a_$(a)_b_$(b)_c_$(c).jld2"
    if isfile(fileName_c)
        println("loading "*fileName_c)
        c_iipDyn_mat_subst = 1.0 .* load_object(fileName_c)
    else
        println("substituting c_iipDyn_mat with a=$(a), b=$(b), c=$(c)")
        fileName_c2 = "CaseStudy/$(lattice_type)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_J1_$(1*j1)_J2_$(1*j2)_J3_$(1*j3)_J4_$(1*j4).jld2"
        c_iipDyn_mat = load_object(fileName_c2)
        c_iipDyn_mat_subst = 1.0 .* get_c_iipDyn_mat_subst(c_iipDyn_mat,hte_lattice,a,b,c);
        save_object(fileName_c,c_iipDyn_mat_subst)
    end

    return c_iipDyn_mat_subst
end

function get_QFI(moments)
    if length(moments) < 7
        moments = vcat(moments, zeros(Int, 7 - length(moments)))
    end
    x_poly = Polynomial([0, 1]) 
    f_Q_poly = 
    x_poly^2    * moments[2] -
    x_poly^4/12 * moments[3] +
    x_poly^6/120 * moments[4] -
    (34 * x_poly^8 / 40320) * moments[5] +
    (62 * x_poly^10 / 725760) * moments[6] -
    (1382 * x_poly^12 / 159667200) * moments[7]
    return f_Q_poly
end


if false
    #calculate the moments
    a = 0.12
    sc = sqrt(1+a^2)
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a)
    k = (4*pi/3,0)
    c_kDyn_mat = get_c_k([k],c_iipDyn_mat_subst,hte_lattice)[1]
    m_vec = get_moments_from_c_kDyn(c_kDyn_mat)#[1:7]
    #println(length(m_vec))
    f_Q_poly = get_QFI(m_vec)

    x0_vec = collect(0.1:0.01:8.0) #J/T
    x_vec = x0_vec./sc # J1/T

    f = 0.25 #0.25
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
    u_vec = tanh.(f .* x_vec)
    coeffs_x = coeffs(f_Q_poly)
    p_u = Polynomial(ufromx_mat*coeffs_x)

    plt = plot(xlabel="J/T", ylabel=L"f_Q[k,T]", title="$(lattice_type), a=$(a), k=$(k)", legend=:bottomright)
    Plots.plot!(plt,x0_vec,robustpade(p_u,5,6).(u_vec),color="blue",linestyle=:dashdot,alpha=0.7,label="u-Padé[5,6] (f=$f)")
    Plots.plot!(plt,x0_vec,robustpade(p_u,4,5).(u_vec),color="red",linestyle=:dash,alpha=0.7,label="u-Padé[4,5] (f=$f)")
    
    #display(plt)
    #savefig(plt,"Images/QFI/$(lattice_type)_a_$(a).png")

    inv_x0_vec = 1.0 ./ x0_vec
    inv_x_vec = 1.0 ./ x_vec
    
    plt2 = plot(xlabel="T/J", ylabel=L"f_Q[k,T]", title="J2/J1=$(a)", xlims=(0,1.25), ylims=(0,5.05), xlabelfontsize = 14, xtickfontsize = 12, ylabelfontsize = 14, ytickfontsize = 12, legendfontsize = 12, titlefontsize=14)

    # comparison data
    if lattice_type == "triang"
        x_sh = collect(0.1:0.1:2.0)./sc
        if a == 0.0
            y_sh = [3.706,2.397,1.667,1.314,1.079,0.929,0.778,0.683,0.595,0.536,0.485,0.440,0.397,0.357,0.324,0.298,0.273,0.252,0.237,0.221] # N=30
            y_sh2 = [3.544,2.470,1.746,1.346,1.089,0.922,0.784,0.683,0.595,0.536,0.485,0.440,0.397,0.357,0.324,0.298,0.273,0.252,0.237,0.221] # N=24
            y_sh3 = [3.192,2.451,1.759,1.384,1.149,0.966,0.829,0.717,0.630,0.556,0.486,0.440,0.397,0.357,0.324,0.298,0.273,0.252,0.237,0.221] # N=18
            col = "orange"

            Plots.plot!(plt2,[NaN],[NaN],color=:white,label=" ",linealpha=0,marker=:none)
            Plots.plot!(plt2,[NaN],[NaN],color=:white,label="[Shimokawa...2025]",linealpha=0,marker=:none)
            Plots.plot!(plt2,x_sh,y_sh,color=col,marker=:circle,markersize=6,alpha=0.3,label="N=30")
            Plots.plot!(plt2,x_sh,y_sh2,color=col,marker=:rect,markersize=6,alpha=0.3,label="N=24")
            Plots.plot!(plt2,x_sh,y_sh3,color=col,marker=:utriangle,markersize=6,alpha=0.3,label="N=18")
        
        elseif a == 0.12
            y_sh = [2.619,1.730,1.351,1.110,0.943,0.808,0.713,0.622,0.547,0.486,0.434,0.390,0.353,0.322,0.296,0.269,0.250,0.229,0.214,0.198] # N=30
            y_sh2 = [2.544,1.791,1.367,1.110,0.943,0.808,0.713,0.622,0.547,0.486,0.434,0.390,0.353,0.322,0.296,0.269,0.250,0.229,0.214,0.198] # N=24
            y_sh3 = [2.147,1.830,1.423,1.175,0.992,0.847,0.732,0.641,0.547,0.486,0.434,0.390,0.353,0.322,0.296,0.269,0.250,0.229,0.214,0.198] # N=18
            col = "orange"

            Plots.plot!(plt2,[NaN],[NaN],color=:white,label=" ",linealpha=0,marker=:none)
            Plots.plot!(plt2,[NaN],[NaN],color=:white,label="[Shimokawa...2025]",linealpha=0,marker=:none)
            Plots.plot!(plt2,x_sh,y_sh,color=col,marker=:circle,markersize=6,alpha=0.3,label="N=30")
            Plots.plot!(plt2,x_sh,y_sh2,color=col,marker=:rect,markersize=6,alpha=0.3,label="N=24")
            Plots.plot!(plt2,x_sh,y_sh3,color=col,marker=:utriangle,markersize=6,alpha=0.3,label="N=18")    
        end

    elseif lattice_type == "chain"
        x_kcuf = [6,50,75,150,200,300] .* (0.08617/33.5) #T/J1
        sc2 = (8*spin_length^2)
        y_kcuf = [1.86,1.24,0.96,0.58,0.45,0.20] .* sc2
        x_bethe = [6*(0.08617/33.5)]
        y_bethe = [2.291*sc2]
        Plots.scatter!(plt2,x_kcuf,y_kcuf,color="orange",marker=:rect,markersize=6,alpha=0.3,label=L"KCuF_3") 
        cap = 0.01
        for (xi, yi, ei) in zip(x_kcuf./sc, y_kcuf, [0.14,0.10,0.14,0.09,0.06,0.15] .* (0.5*sc2))
            Plots.plot!([xi, xi], [yi - ei, yi + ei], color=:orange, linewidth=0.8, label=false)
            Plots.plot!([xi - cap, xi + cap], [yi + ei, yi + ei], color=:orange, linewidth=0.8, label=false)
            Plots.plot!([xi - cap, xi + cap], [yi - ei, yi - ei], color=:orange, linewidth=0.8, label=false)
        end
        Plots.scatter!(plt2,x_bethe./sc,y_bethe,color="green",marker=:circle,markersize=6,alpha=0.3,label="Bethe ansatz")
    end
    
    Plots.plot!(plt2,inv_x0_vec,f_Q_poly.(x_vec),color="gray",label="x-Series")
    Plots.plot!(plt2,inv_x0_vec,robustpade(p_u,5,6).(u_vec),color="blue",linestyle=:dash,alpha=0.7,label="u-Padé[5,6] (f=$f)")
    Plots.plot!(plt2,inv_x0_vec,robustpade(p_u,4,5).(u_vec),color="red",linestyle=:dashdot,alpha=0.7,label="u-Padé[4,5] (f=$f)")
    display(plt2)
    #savefig(plt2,"Images/$(lattice_type)_a_$(a).png")

    if false
        plt_sh4 = plot(xlabel="T/J", ylabel=L"f_Q[\mathbf{k},T]", title="$(lattice_type), J2/J1=$(a)", xscale=:log10,yscale=:log10,xlims=(0.08,1), ylims=(0.5,10), xlabelfontsize = 14, xtickfontsize = 12, ylabelfontsize = 14, ytickfontsize = 12, legendfontsize = 11, titlefontsize=14)

        #Plots.plot!(plt_sh4, inv_x_vec, f_Q_poly.(x_vec),color="gray",label="x-Series")
        Plots.plot!(plt_sh4,inv_x0_vec,robustpade(p_u,5,6).(u_vec),color="blue",linestyle=:dash,alpha=0.7,label="u-Padé[5,6] (f=$f)")
        Plots.plot!(plt_sh4,inv_x0_vec,robustpade(p_u,4,5).(u_vec),color="red",linestyle=:dashdot,alpha=0.7,label="u-Padé[4,5] (f=$f)")

        # Powerlaw
        # Model: y = A * x^α
        powerlaw(x, p) = p[1] .* x .^ p[2]  # p[1] = A, p[2] = α

        xdata = collect(1:0.01:2.5) # J1/T
        ydata = robustpade(p_u,4,5).(tanh.(f .* xdata))

        p0 = [1.0, -1.0]  # initial guess [A, α]
        fit = curve_fit(powerlaw, xdata, ydata, p0)

        A, α = fit.param
        println("Fitted A = $A")
        println("Fitted alpha = $α")
        # Generate smooth x values for fit line
        xfine = range(0.4, stop=10, length=200) # J1/T
        yfit = powerlaw(xfine, fit.param)

        # Plot fitted curve
        Plots.plot!((xfine.*sc).^(-1), yfit, label="Fit: y = $(round(A, digits=2))·x^($(round(α, digits=2)))", lw=2, color=:green)
        
        display(plt_sh4)
        #savefig(plt_sh4,"Images/$(lattice_type)_$(a)_log.png")
    end
end

# a fixed, vary k
if false
    a = 1.1
    sc = sqrt(1+a^2)
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a)
    
    plt3 = plot(xlabel="T/J", ylabel=L"f_Q[k,T]", title="J2/J1=$(a)", xlims=(0,2.05), ylims=(0,4), xlabelfontsize = 14, xtickfontsize = 12, ylabelfontsize = 14, ytickfontsize = 12, legendfontsize = 12, titlefontsize=14) #,xlims=(0,2.05), ylims=(0,4)
    
    x0_vec = collect(0.1:0.01:7.0) #J/T
    x_vec = x0_vec./sc # J1/T
    inv_x0_vec = 1.0 ./ x0_vec #T/J
    inv_x_vec = 1.0 ./ x_vec #T/J1

    f_k = [0.4,0.3,0.5]
    m_k = [6,5,5]
    for (i,k) in enumerate([(pi,pi),(pi,0),(pi/2,pi/2)]) #[(pi,pi),(pi,0),(pi/2,pi/2)]
        #f = f_k[i]
        #m = m_k[i]
        f = 0.4
        m = 5
        ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
        u_vec = tanh.(f .* x_vec)
        c_kDyn_mat = get_c_k([k],c_iipDyn_mat_subst,hte_lattice)[1]
        m_vec = get_moments_from_c_kDyn(c_kDyn_mat)
        f_Q_poly = get_QFI(m_vec)
        coeffs_x = coeffs(f_Q_poly)
        p_u = Polynomial(ufromx_mat*coeffs_x)
        #Plots.plot!(plt3, inv_x0_vec, f_Q_poly.(x_vec),color="gray",label="x-Series")
        Plots.plot!(plt3,inv_x0_vec,robustpade(p_u,m,6).(u_vec),linestyle=:dash,alpha=0.7,label="k=($(round(k[1],digits=2)),$(round(k[2],digits=2))) , u-Padé[5,6] (f=$f)")
        Plots.plot!(plt3,inv_x0_vec,robustpade(p_u,4,5).(u_vec),linestyle=:dashdot,alpha=0.7,label="k=($(round(k[1],digits=2)),$(round(k[2],digits=2))), u-Padé[4,5] (f=$f)")
    end
    
    display(plt3)
    #savefig(plt3,"Images/$(lattice_type)_a_$(a).png")
end

# vary a, fixed k
if false
    f_k = [0.4,0.2,0.3,0.3,0.4]
    m_k = [5,5,5,5,6]
    k = (4*pi/3,0) #(pi,pi) #(pi/2,pi/2)

    plt4 = plot(xlabel="T/J", ylabel=L"f_Q[k,T]", title="k=(π,π)", xlims=(0,2.05), ylims=(0,4), xlabelfontsize = 14, xtickfontsize = 12, ylabelfontsize = 14, ytickfontsize = 12, legendfontsize = 12, titlefontsize=14) #,xlims=(0,2.05), ylims=(0,4)
    
    for (i,a) in enumerate([0.0,0.047,0.08,0.12,0.25,0.5,1.0,2.0]) #[0.2,1.1,1.25,1.4,4.0] [0.0,0.047,0.08,0.12,0.25,0.5,1.0,2.0]
        sc = sqrt(1+a^2)
        x0_vec = collect(0.1:0.01:4.0) #J/T
        x_vec = x0_vec./sc # J1/T
        inv_x0_vec = 1.0 ./ x0_vec #T/J
        inv_x_vec = 1.0 ./ x_vec # T/J1
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a)
        #f = f_k[i]
        #m = m_k[i]
        f = 0.23
        m = 5
        ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
        u_vec = tanh.(f .* x_vec)
        c_kDyn_mat = get_c_k([k],c_iipDyn_mat_subst,hte_lattice)[1]
        m_vec = get_moments_from_c_kDyn(c_kDyn_mat)
        f_Q_poly = get_QFI(m_vec)
        coeffs_x = coeffs(f_Q_poly)
        p_u = Polynomial(ufromx_mat*coeffs_x)
        #Plots.plot!(plt4, inv_x0_vec, f_Q_poly.(x_vec),color="gray",label="x-Series")
        Plots.plot!(plt4,inv_x0_vec,robustpade(p_u,m,6).(u_vec),linestyle=:dash,alpha=0.7,label="a=$(a)") # , u-Padé[5,6] (f=$f)")
        #Plots.plot!(plt4,inv_x0_vec,robustpade(p_u,4,5).(u_vec),linestyle=:dashdot,alpha=0.7,label="a=$(a)", u-Padé[4,5] (f=$f)")
    end
    
    display(plt4)
    #savefig(plt4,"Images/$(lattice_type)_k_$(k).png")
end

# BZ-plot of nQFI
if false
    a = 2.0
    sc = sqrt(1+a^2)
    x0 = 2.0 # J/T^
    x = x0/sc # J1/T
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a)
    f = 0.7
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
    u = tanh.(f .* x)
    pade_order=[5,6]

    LL=60
    step_vec = collect(-3*LL:1:3*LL)/LL
    function kk(f1,f2)
        #f1 .* G1 + f2 .* G2
        f1 .* (π,0) .+ f2 .* (0,π) 
    end
    k_vec = [kk(f1,f2) for f1 in step_vec for f2 in step_vec if abs(kk(f1,f2)[1]/π)<=1.5 && abs(kk(f1,f2)[2]/π)<=1.5]
    kx_vec, ky_vec = [k[1] for k in k_vec],[k[2] for k in k_vec]

    coeffs_x_vec = [coeffs(get_QFI(get_moments_from_c_kDyn(get_c_k(k, c_iipDyn_mat_subst,hte_lattice)))) for k in k_vec]

    if f==0.0 ### use x-Polynomial
        z_data = [get_pade(Polynomial(coeffs_x),pade_order...)(x) for coeffs_x in coeffs_x_vec]
    else ### use u-Polynomial
        z_data = [get_pade(Polynomial(ufromx_mat*coeffs_x),pade_order...)(u) for coeffs_x in coeffs_x_vec]
        #z_data = [Polynomial(ufromx_mat*coeffs_x)(u) for coeffs_x in coeffs_x_vec]
    end
end

if false
    ### prepare plot
    using CairoMakie
    xPlots,yPlots=1,1
    fig = Figure(size=(1*aps_width*xPlots,1*0.8*aps_width*yPlots),fontsize=8)
    ax=Axis(fig[1,1],xlabel=L"k_x/\pi",ylabel=L"k_y/\pi",xlabelsize=9,ylabelsize=9,aspect=1,title="$(lattice_type): Padé$pade_order x="*string(x0)* " J2/J1=$a f=$f")
    hm=CairoMakie.heatmap!(ax,kx_vec/π,ky_vec/π,z_data,colormap=:viridis,colorrange=(0.0,1.0),highclip=:white,lowclip=:black)
    lines!(ax,[-1,-1,1,1,-1],[-1,1,1,-1,-1],color=:grey)
    #lines!(ax,[4/3*cos(α*π/3) for α in 0:6],[4/3*sin(α*π/3) for α in 0:6],color=:grey)
    CairoMakie.Colorbar(fig[:, end+1], hm, label="nQFI",)
    CairoMakie.ylims!(ax,(-1.5,1.5))
    CairoMakie.xlims!(ax,(-1.5,1.5))
    resize_to_layout!(fig)
    display(fig)
    save("Images/"*string(lattice_type)*"_a$a"*"_x$(x0).png",fig)
end
