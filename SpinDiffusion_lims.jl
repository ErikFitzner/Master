using JLD2, DelimitedFiles, SimpleWeightedGraphs, Plots, Richardson, CairoMakie, Trapz, Statistics, QuadGK

include("plotConventions.jl")
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 


### load graph evaluations
spin_length = 1/2
bc1 = spin_length*(spin_length+1)/3
n_max = 14

### prepare lattice
lattice_type = "aniso_square" #"Shastry-Sutherland", "ladder_collapse"
d = 2
ex = zeros(d)  # creates a d-dimensional zero vector
ex[1] = 1      # set first component to 1

# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = true
j3 = false
j4 = false

L = 14

hte_lattice = getLattice(L,lattice_type,j1,j2,j3,j4);  #"chain", "square"

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

""" For translation-invariant, sufficiently isotropic and inversionsymmetric systems, the k→0 limit can be obtained directly from the Dyn-HTE expansion of Gii′(iνm).
The way we take the k→0 limit sets the transport direction, here: x-direction """
function get_moments_Gklim(c_iipDyn_mat, hte_lattice, ex)
    lattice = hte_lattice.lattice
    center_sites = hte_lattice.basis_positions
    z = zeros(size(c_iipDyn_mat[1]))
    for b in 1:length(center_sites)
        for i in 1:size(c_iipDyn_mat)[1]
            z += (dot(ex, getSitePosition(lattice,i) .- getSitePosition(lattice,center_sites[b])))^2 * c_iipDyn_mat[i,b]
        end
    end
    Gklim = (-1/2) * z / length(center_sites)
    moments = get_moments_from_c_kDyn(Gklim)
    return moments[2:end]
end

function Diffusion_Gauß(Tχ, moments,x)
    m2 = moments[1](x)
    m4 = moments[2](x)
    result = sqrt(pi)/(Tχ(x)*sqrt(2))*sqrt(m2^3/m4)
    return result
end

""" get σ(w) from δ_vec at w=ω/J and broadening η"""
function sigma(δ_vec::Vector{Float64},x::Float64,w::Float64,η::Float64)::Float64
    #res = 1/pi * real(contFrac(1im * w + η,δ_vec))
    res = 1/pi * real(contFrac(1im * w + η,δ_vec.^2))
    if x==0.0
        return pi * res
    else
        return pi * x * res
    end
end

""" get σ(w) from δ_vec at w=ω/J and broadening η"""
function sigma_neu(δ_vec::Vector{Float64},x::Float64,w::Float64,extrap_params::Vector{Float64})::Float64
    a = extrap_params[1]
    b = extrap_params[2] + (length(δ_vec)-1)*a 
    if a < 0.0
        println("Slope of δ extrapolation must be non-negative")
        #throw(ErrorException("Slope of δ extrapolation must be non-negative"))
        a = abs(a)
    end

    res = 1/pi * real(contFracwithTerminator(1im * w,δ_vec,a,b))#^2
    if x==0.0
        return pi * res
    else
        return pi * x * res
    end
end

function get_upade(Tχ::Polynomial{Float64, :x};f::Float64=0.7,m::Int=6,n::Int=5)
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
    coeffs_x = coeffs(Tχ)
    p_u = Polynomial(ufromx_mat*coeffs_x)
    return get_pade(p_u,m,n)
end

###################################################################
###### Analytic continuation ######################################
###################################################################

function get_delta_vec_ext(x0::Float64,moments;f::Float64=0.7,r_max::Int64=3)
    n_max_d = n_max-2
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max_d+1,f)
    poly_x = Polynomial([0,1],:x)
    u0 = tanh.(f .* x0)
    if x0 == 0
        m0 = []
        for m_idx=1:r_max+1
            push!(m0, get_pade(moments[m_idx],7-m_idx,6-m_idx))
        end
        δ_vec,r_vec = fromMomentsToδ([m(x0) for m in m0])
    else
        m0 = Float64[]
        for r in 0:r_max
            xm_norm_r = coeffs(poly_x * (moments[1+r]/moments[1+r](0)))
            p_u = Polynomial(ufromx_mat[1:n_max_d+2-2*r,1:n_max_d+2-2*r]*xm_norm_r)
            push!(m0,moments[1+r](0)/x0 * get_pade(p_u,7-r,6-r)(u0))
        end
        δ_vec,r_vec = fromMomentsToδ(m0)
    end

    δ_vec_ext = extrapolate_δvec(δ_vec,length(δ_vec)-2,length(δ_vec)-1,2000,false)
    b_vec = [sqrt(d) for d in δ_vec]
    b_vec_ext = extrapolate_δvec(b_vec,length(b_vec)-2,length(b_vec)-1,2000,false)

    #println(b_vec)

    ##### check moments #####
    if false
        x_vec = 0.01:0.01:2.5
        m_vec = moments
        m_vec_times_x = [m_vec[i]*Polynomial([0,1]) for i=1:length(m_vec)]

        ##EXTRAPLOATION OF MOMENTS
        #basic pade 
        m_vec_extrapolated_pade = []
        for m_idx=1:length(m_vec)-2
            push!(m_vec_extrapolated_pade, extrapolate_series(m_vec[m_idx],"pade",(7-m_idx,7-m_idx)))
        end
        m_vec_times_x_extrapolated_u_pade1 = []
        m_vec_times_x_extrapolated_u_pade2 = []
        for m_idx=1:length(m_vec)-2
            push!(m_vec_times_x_extrapolated_u_pade1, extrapolate_series(m_vec_times_x[m_idx],"u_pade",(8-m_idx,7-m_idx,f)))
            push!(m_vec_times_x_extrapolated_u_pade2, extrapolate_series(m_vec_times_x[m_idx],"u_pade",(7-m_idx,6-m_idx,f)))
        end

        #plot the moments 
        plt_m = Plots.plot([0],[0],label="",xlabel=L"x",ylabel=L"x \cdot m_{d,2r}(x)/m_{d,2r}(0)",legend=:topleft,xlim=(-0.2,2.5),size=(1.5*aps_width,aps_width))
        Plots.plot!(plt_m,[0],[0],label="x bare",color = "grey",linestyle = linestyle_vec[1],linewidth=0.4)
        Plots.plot!(plt_m,[0],[0],label="u Padé [8-r,7-r]",color = "grey",linestyle = linestyle_vec[2],alpha=0.5)
        Plots.plot!(plt_m,[0],[0],label="u Padé [7-r,6-r]",color = "grey",linestyle = linestyle_vec[3])
        #Plots.plot!(plt_m,[0],[0],label="Padé [7-r,7-r]",color = "grey",linestyle = linestyle_vec[4])
        for i=1:4
            Plots.plot!(plt_m,[0],[0],label="r="*string(i),color = color_vec[i])
            Plots.plot!(plt_m,x_vec[1:180],m_vec_times_x[i].(x_vec[1:180])./(m_vec[i](0)),label = "",alpha= 0.7,color = color_vec[i],linestyle = linestyle_vec[1],linewidth=0.5)
            Plots.plot!(plt_m,x_vec,m_vec_times_x_extrapolated_u_pade1[i].(tanh.(f.*x_vec))./(m_vec[i](0)),label = "",alpha= 0.5,color = color_vec[i],linestyle = linestyle_vec[2],linewidth=0.5)
            Plots.plot!(plt_m,x_vec,m_vec_times_x_extrapolated_u_pade2[i].(tanh.(f.*x_vec))./(m_vec[i](0)),label = "",alpha= 1,color = color_vec[i],linestyle = linestyle_vec[3],linewidth=0.5)
            #Plots.plot!(plt_m,x_vec,x_vec.*m_vec_extrapolated_pade[i].(x_vec)./(m_vec_extrapolated_pade[i](0.0001)),label = "",alpha= 1,color = color_vec[i],linestyle = linestyle_vec[4],linewidth=0.5)
        end
        #savefig(plt_m,"Images/moments_extrapolation_$(lattice_type).svg")
        #SHOW PLOT
        ##############
        display(plt_m)
    end

    ##### check extrapolation #####
    if false
        plt_δ = Plots.plot([0],[0],label="",xlabel=L"r_\text{max}",legend=:topleft,size=(1.5*0.8*aps_width,0.8*aps_width))
        Plots.scatter!(plt_δ,r_vec,δ_vec,label="",color=thermalCol13_vec[1])
        Plots.scatter!(plt_δ,r_vec,b_vec,label="",color=thermalCol13_vec[13])
        Plots.plot!(plt_δ,0:7,δ_vec_ext[0+1:8],label="δ_vec",color=thermalCol13_vec[1])
        Plots.plot!(plt_δ,0:7,b_vec_ext[0+1:8],label="b_vec",color=thermalCol13_vec[13])
        Plots.xticks!(0:7,string.(1:8))
        display(plt_δ)
        #savefig(plt_δ,"Images/bVec_$(lattice_type)_x0_$(x0)_rmax_$(r_max).png")
    end
    return δ_vec_ext, b_vec_ext
end

function get_delta_vec_ext_neu(x0::Float64,moments;f::Float64=0.7,r_max::Int64=3,r_min::Int64=3)
    intercept0 = true
    n_max_d = n_max - 2
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max_d+1,f)
    poly_x = Polynomial([0,1],:x)
    u0 = tanh.(f .* x0)
    #m0 = [moments[i](x0) for i in 1:length(moments)]
    if true
        m0 = Float64[]
        for r in 0:r_max
            xm_norm_r = coeffs(poly_x * (moments[1+r]/moments[1+r](0)))
            p_u = Polynomial(ufromx_mat[1:n_max_d+2-2*r,1:n_max_d+2-2*r]*xm_norm_r)
            push!(m0,moments[1+r](0)/x0 * get_pade(p_u,6-r,5-r)(u0))
        end
    end
    δ_vec,r_vec = fromMomentsToδ(m0)
    
    δ_vec = [sqrt(d) for d in δ_vec]
    
    # find last index where δ_vec is non-negative
    idx = findfirst(<(0), δ_vec)
    lastidx = isnothing(idx) ? length(δ_vec) : idx - 1
    r_max_eff = min(lastidx - 1, r_max)
    #r_min_eff = min(r_min, r_max_eff)
    #r_min_eff = max(0, r_max_eff-1)
    r_min_eff = r_max_eff

    if r_max_eff < 1
        println("WARNING: negative δ1, putting δ0 = 0")
        extrap_params = [1.0,0.0]
        δ_vec = [0.0,1.0]
        r_max_eff = 0
    else
        #println("r_max_eff = $(r_max_eff), r_min_eff = $(r_min_eff)")
        extrap_params = get_extrapolation_params(δ_vec[1:r_max_eff+1],r_min_eff,r_max_eff,intercept0)
    end
    
    #println(extrap_params)
    #println(δ_vec)
    
    # plot delta
    if false
        plt_δ = Plots.plot([0],[0],title="$(lattice_type)",label="",xlabel=L"r",ylabel=L"b_r",legend=:bottomright)
        Plots.scatter!(plt_δ,r_vec,δ_vec,label="")
        δ_vec_ext = extrap_params[2] .+ extrap_params[1] .* (0:8)
        Plots.plot!(plt_δ,0:7,δ_vec_ext[0+1:8],label="",color=thermalCol13_vec[1])  # r_max+1
        display(plt_δ)
        #savefig(plt_δ,"Images/bVec_$(lattice_type)_x0_$(x0)_rmax_$(r_max).png")
    end

    return δ_vec, extrap_params
end


"""
    fermionic_integral(beta, J)

Evaluates ∫_{-π}^{π} dk sin²(k) f(ε_k)(1 - f(ε_k)),
with ε_k = J cos k and f(ε) = 1 / (exp(β ε) + 1).
"""
function fermionic_integral(x::Float64)::Float64
    integrand(k) = begin
        f = 1.0 / (exp(x * cos(k)) + 1.0)
        sin(k)^2 * f * (1.0 - f)
    end
    val, err = quadgk(integrand, -pi, pi, rtol=1e-10)
    return val
end

function Tχ_exact(x::Float64)::Float64
    integrand(k) = begin
        f = 1.0 / (exp(x * cos(k)) + 1.0)
        f * (1.0 - f)
    end
    val, err = quadgk(integrand, -pi, pi, rtol=1e-10)
    return 1/(2*pi) * val
end


###################################################################
# comparison σ(w) for different values of x at fixed a
function sigma_w_x(x_vec::Vector{Float64},w_vec::Vector{Float64},a::Float64;r_max::Int64=3,f::Float64=0.7,b::Float64=0.0,c::Float64=0.0)
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
    moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
    sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS
    
    plt = Plots.plot([0],[0],label="",xlabel=L"w",ylabel=L"\sigma(w)",title="Spin conductivity, $(lattice_type), J2/J1=$(a)",size=(1.5*aps_width,aps_width), legend=:topright)

    for x in x_vec # J/T
        x0 = x/sc  # J1/T

        δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r_max,f=f)[2]
        sigma_ω = [sigma(δ_vec_ext,x0,w*sc,0.001) for w in w_vec]

        #δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)
        #sigma_ω = [sigma_neu(δ_vec,x0,w*sc,extrap_params) for w in w_vec]

        Plots.plot!(plt, w_vec, sigma_ω, label="x=$(x)")
    end
    #display(plt)
    #savefig(plt,"Images/SpinConductivity_$(lattice_type)_a_$(a).png")
end

# comparison D_S/J(w) for different values of x at fixed a
function DS_w_x(x_vec::Vector{Float64},w_vec::Vector{Float64},a::Float64;r_max::Int64=3,f::Float64=0.7,b::Float64=0.0,c::Float64=0.0)
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
    moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
    c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
    Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
    Tχ_upade = get_upade(Tχ,f=f)

    sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS
    
    plt = Plots.plot([0],[0],label="",xlabel=L"w",ylabel=L"D_S/J(w)",size=(0.8*1.5*aps_width,0.8*aps_width), legend=:topright) # title="$(lattice_type) lattice, J2/J1=$(a)",

    for x in x_vec # J/T
        x0 = x/sc  # J1/T
        
        for r in 1:r_max # 2_r_max
            δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r,f=f)[2]
            
            if x0 == 0.0
                DS_w = [1/Tχ_upade(tanh(f*x0)) * sigma(δ_vec_ext,x0,w*sc,0.001) for w in w_vec]
            else
                DS_w = [1/x0 * 1/Tχ_upade(tanh(f*x0)) * sigma(δ_vec_ext,x0,w*sc,0.001) for w in w_vec]
            end

            Plots.plot!(plt, w_vec, DS_w, label=L"r_{\text{max}}=" * "$(r+1)", color=color_vec[r]) # , label="x=$(x)"
        end
        
    end
    display(plt)
    #savefig(plt,"Images/DS_$(lattice_type)_x_$(x_vec[1]).svg")  # _a_$(a)
end

# comparison comparison σ(w) for different values of a at fixed x
function sigma_w_a(x::Float64, w_vec::Vector{Float64}, a_vec::Vector{Float64};b::Float64=0.0,c::Float64=0.0,r_max::Int64=3,f::Float64=0.7)
    plt2 = Plots.plot([1],[1],label="",xlabel=L"w=\omega/J_x",ylabel=L"D_{s,x}(w)/J_x",size=(1.2*aps_width,0.8*aps_width), legend=:bottomleft, xscale=:log10, yscale=:log10) # ,title="Spin conductivity, $(lattice_type), x=$(x)"
    ϵ = 1e-8
    sigma_w_0 = []
    for a in a_vec
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
        moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
        c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
        Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
        Tχ_upade = get_upade(Tχ,f=f)
        sc = 1 #sqrt(1+a^2+b^2+c^2) # scale to w/J and JS
        x0 = x/sc # J1/T
        
        δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r_max,f=f)[2]
        sigma_w = [sigma(δ_vec_ext,x0,w*sc,ϵ) for w in w_vec]

        if x0 == 0.0
            DS_w = [1/Tχ_upade(tanh(f*x0)) * sigma(δ_vec_ext,x0,w*sc,ϵ) for w in w_vec]
        else
            DS_w = [1/x0 * 1/Tχ_upade(tanh(f*x0)) * s for s in sigma_w]
        end

        Plots.plot!(plt2, w_vec, DS_w, label=L"J_y/J_x=" * "$(a)")  # sigma_w
        push!(sigma_w_0, DS_w[1])
    end

    if lattice_type == "aniso_square"
        sigma_exact = [x/2 * fermionic_integral(x) * (1/pi * ϵ/(w^2+ϵ^2)) for w in w_vec]
        DS_exact = [1/2 * fermionic_integral(x) * 1/Tχ_exact(x) * (1/pi * ϵ/(w^2+ϵ^2)) for w in w_vec]
        Plots.plot!(plt2, w_vec, DS_exact, linestyle=:dash, color=:black, label="XX chain exact")
    end
    display(plt2)
    #savefig(plt2,"Images/SpinConductivity_$(lattice_type)_x_$(x).svg")

    if false
        # --- log–log fit ---
        mask = (a_vec .< 0.11) .& (sigma_w_0 .> 0)
        x_fit = log.(a_vec[mask])
        y_fit = log.(sigma_w_0[mask])
        m = cov(x_fit, y_fit) / var(x_fit)
        b = mean(y_fit) - m * mean(x_fit)
        α = -m
        C = exp(b)
        x_vec = collect(minimum(a_vec):0.001:maximum(a_vec))
        fit_curve = C .* x_vec.^(-α)
        #println(C)

        plt3 = Plots.plot([1],[1],label="",xlabel=L"J_y/J_x",ylabel=L"D_{s,x}(w\to 0)/J_x",size=(1.2*aps_width,0.8*aps_width), legend=:topright, xscale=:log10, yscale=:log10) # , xlims=(0,0.21)
        Plots.scatter!(plt3, a_vec, sigma_w_0, marker=:circle, label="Dyn-HTE")
        Plots.plot!(plt3, x_vec, fit_curve, linestyle=:dash, color=:black, label = L"y = %$(round(C, sigdigits=2))\, (J_y/J_x)^{-%$(round(α,digits=2))}")
        display(plt3)
        #savefig(plt3,"Images/SpinConductivity_$(lattice_type)_x_$(x)_w_0.png")
    end
end

function DS_w_a_moments(x::Float64, a_vec::Vector{Float64};b::Float64=0.0,c::Float64=0.0)
    plt2 = Plots.plot([1],[1],label="",size=(aps_width,0.8*aps_width),xlabel=L"r",ylabel=L"m_{d,2r}"*"(x=$(x))",legend=:topright, yscale=:log10)
    for a in a_vec
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
        moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
        #println(moments[1](0.0))
        Plots.scatter!(plt2, 1:length(moments), [moments[i](x) for i in 1:length(moments)], label=L"J_y/J_x="*"$a")
    end
    display(plt2)
    #savefig(plt2,"Images/DS_$(lattice_type)_moments.png")
end


function DS_w_a_int(x::Float64, w_vec::Vector{Float64}, a_vec::Vector{Float64};b::Float64=0.0,c::Float64=0.0,r_max::Int64=3,f::Float64=0.7)
    plt2 = Plots.plot([0],[0],label="",xlabel=L"w_c",ylabel=L"\int_{-w_c}^{w_c} D_{s,x}(w) dw/J_x",size=(0.67*aps_width,0.5*aps_width), legend=:bottomright, yticks=([0.0,pi/4,pi/2],["0.0","π/4","π/2"]),title="x=$x") #, ylims=(-0.01,0.29))
    for a in a_vec
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
        moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
        c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
        Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
        Tχ_upade = get_upade(Tχ,f=f)
        sc = 1 #sqrt(1+a^2+b^2+c^2) # scale to w/J and JS
        x0 = x/sc # J1/T

        δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r_max,f=f)[2]

        if x0 == 0.0
            DS_w = [1/Tχ_upade(tanh(f*x0)) * sigma(δ_vec_ext,x0,w*sc,1e-8) for w in w_vec]
        else
            DS_w = [1/x0 * 1/Tχ_upade(tanh(f*x0)) * sigma(δ_vec_ext,x0,w*sc,1e-8) for w in w_vec]
        end

        wc_vec = abs.(w_vec[w_vec .>= 0])
        W = zeros(length(wc_vec))

        for (i, wc) in enumerate(wc_vec)
            mask = abs.(w_vec) .<= wc
            W[i] = 2 * trapz(w_vec[mask], DS_w[mask])
        end
        Plots.plot!(plt2, wc_vec, W, label=L"J_y/J_x="*"$a")
    end
    println(Tχ_exact(x))
    Plots.hline!(plt2, [fermionic_integral(x)/(2*Tχ_exact(x))], linestyle=:dash, color=:black, label="XX chain exact")
    #Plots.hline!(plt2, [1/4], linestyle=:dash, color=:black, label="Exact")
    display(plt2)
    #savefig(plt2,"Images/DS_$(lattice_type)_r_x_$(x)_int.svg")
end

# comparison D_S with Gauß and literature
function comp_DS_Dyn_Gauß(x_vec::Vector{Float64},a::Float64,str::String;b::Float64=0.0,c::Float64=0.0,r_max=Int64=3,f::Float64=0.7)
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
    moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
    c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
    Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
    Tχ_upade = get_upade(Tχ,f=f)
    sc = 1 # sqrt(1+a^2+b^2+c^2) # scale to w/J and JS

    y_vec = [Diffusion_Gauß(Tχ, moments, x/sc) for x in x_vec]
    plt3 = Plots.plot([0],[0],label="",xlabel=L"J/T",ylabel=L"D_S/J",title="Spin Diffusion, $(lattice_type), J2/J1=$(a)",size=(1.5*aps_width,aps_width), legend=:bottomleft)  # title="Spin Diffusion, XX $(lattice_type), " * L"J_\perp/J=1.0"
    Plots.plot!(plt3, x_vec, y_vec, label="Dyn-HTE: Gauß Approximation")

    #sig_x = [1/(i/sc) * 1/Tχ(i/sc) * sigma(get_delta_vec_ext(i/sc,moments,r_max=r_max,f=f),i/sc,w,0.02) for i in x_vec]
    
    sig_x = []
    for x in x_vec #J/T
        x0 = x/sc # J1/T
        δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r_max,f=f)[2]
        
        if x0 == 0.0
            push!(sig_x, 1/Tχ_upade(tanh(f*x0)) * sigma(δ_vec_ext,x0,0.0,0.001))
        else
            push!(sig_x, 1/x0 * 1/Tχ(x0) * sigma(δ_vec_ext,x0,0.0,0.001))
        end
    end
    
    Plots.plot!(plt3, x_vec, sig_x./sc, label="Dyn-HTE")
    #println(sig_x[1]-y_vec[1])

    #### hypercubic lattice
    if str == "hypercubic" && a==0.0
        Plots.scatter!(plt3, [0], [sqrt(pi*bc1/(4*d-2-1/(4*bc1)))], label="[Kopietz1993]", markersize=6)
    elseif str == "square" && a==0.0
    #### square
        Plots.scatter!(plt3, [0], [sqrt(pi*bc1/(4*d-2-1/(4*bc1)))], label="[Kopietz1993]", markersize=6)
        Plots.plot!(plt3,x_vec, sqrt(pi/20) .* (1 .+ 21/20 .* x_vec), linestyle=:dash, label="[SokolPRL1993]")
    else
        println("No literature values")
    end
    display(plt3)
    #savefig(plt3,"Images/SpinDiffusion_$(lattice_type)_lit_a_$(a).png")
end

# comparison D_S(x) for different a
function DS_x_a(x_vec::Vector{Float64},w::Float64,a_vec::Vector{Float64};b::Float64=0.0,c::Float64=0.0,r_max=Int64=3,f::Float64=0.7)
    plt4 = Plots.plot([0],[0],label="",xlabel=L"J_x/T",ylabel=L"D_{s,x}/J_x",size=(0.8*1.5*aps_width,0.8*aps_width), legend=:topleft, xlims=(-0.05,2.05), ylims=(0.65,3.05))  # title="Spin Diffusion, $(lattice_type)",
    for a in a_vec
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
        moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
        c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
        Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
        Tχ_upade = get_upade(Tχ,f=f)
        sc = 1 #sqrt(1+a^2+b^2+c^2) # scale to w/J and JS

        #sig_x = [1/(i/sc) * 1/Tχ(i/sc) * sigma(get_delta_vec_ext(i/sc,moments,r_max=r_max,f=f),i/sc,w,0.02) for i in x_vec]
        
        sig_x = []
        for x in x_vec #J/T
            x0 = x/sc #J1/T

            #δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)
            #push!(sig_x, 1/x0 * 1/Tχ_upade(tanh(f*x0)) * sigma_neu(δ_vec,x0,w*sc,extrap_params))

            δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r_max,f=f)[2] # r_max

            if x0 == 0.0
                push!(sig_x, 1/Tχ_upade(tanh(f*x0)) * sigma(δ_vec_ext,x0,w*sc,0.001))
            else
                push!(sig_x, 1/x0 * 1/Tχ_upade(tanh(f*x0)) * sigma(δ_vec_ext,x0,w*sc,0.001))
            end
        end
        Plots.plot!(plt4, x_vec, sig_x./sc, label=L"J_2/J_1=" * "$(a)")  # label="J2/J1=$(a)"
    end
    
    # Korea
    #=
    #Plots.hline!(plt4, [0.93], linestyle=:dash, color=:gray, label=false)
    Plots.plot!(plt4, [0.65,0.65],[0.5,0.9], linestyle=:dash, color=:gray, label=false)
    Plots.plot!(plt4, [-0.05,0.65],[0.9,0.9], linestyle=:dash, color=:gray, label=false)
    xs, labels = xticks(plt4)[1]
    new_ticks = vcat(xs, 0.65)
    Plots.xticks!(plt4, sort(new_ticks))
    =#
    display(plt4)
    # savefig(plt4,"Images/SpinDiffusion_$(lattice_type).svg")
end

function DS_a_x(x_vec::Vector{Float64},w::Float64,a_vec::Vector{Float64};b::Float64=0.0,c::Float64=0.0,r_max=Int64=3,f::Float64=0.7)
    colors = cgrad(:thermal, length(x_vec), categorical=true)
    plt4 = Plots.plot([1],[1],label="",xlabel=L"J_y/J_x",ylabel=L"D_{s,x}/J_x",size=(aps_width,0.67*aps_width), legend=:topright, xscale=:log10, yscale=:log10) #, ylims=(0.65,1.25))  # title="Spin Diffusion, $(lattice_type)",
    for (i,x) in enumerate(x_vec)
        Ds_list = []
        if x > 0.1
            mask = a_vec .< 2.0
        else
            mask = 1:length(a_vec)
        end
        for a in a_vec[mask]
            c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
            moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
            c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
            Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
            Tχ_upade = get_upade(Tχ,f=f)

            δ_vec_ext = get_delta_vec_ext(x,moments,r_max=r_max,f=f)[2]
            if x == 0.0
                push!(Ds_list, 1/Tχ_upade(tanh(f*x)) * sigma(δ_vec_ext,x,w,1e-8))
            else
                push!(Ds_list, 1/x * 1/Tχ_upade(tanh(f*x)) * sigma(δ_vec_ext,x,w,1e-8))
            end
        end
        Plots.scatter!(plt4, a_vec[mask], Ds_list, color=colors[i], label="J/T=$(round(x,digits=1))")
    end
    fit = 0.5*(3/(2*pi))^0.5 ./ (a_vec .^ 2)  # 0.41
    fit2 = pi^0.5/(2*2^0.5) ./ a_vec  # 0.64
    Plots.plot!(plt4, a_vec, fit, linestyle=:dash, label=L"\frac{1}{2}\sqrt{\frac{3}{2\pi}}\frac{J_x^2}{J_y^2}", color=:black)
    Plots.plot!(plt4, a_vec, fit2, linestyle=:dashdot, label=L"\frac{\sqrt{\pi}}{2\sqrt{2}}\frac{J_x}{J_y}", color=:gray)
    display(plt4)
    savefig(plt4,"Images/SpinDiffusion_$(lattice_type)_r.svg")
end

function DS_x_a_heatmap(x_vec::Vector{Float64},w::Float64,a_vec::Vector{Float64};b::Float64=0.0,c::Float64=0.0,r_max=Int64=3,f::Float64=0.7)
    Z = zeros(length(a_vec),length(x_vec))
    for (i,a) in enumerate(a_vec)
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
        moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
        c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
        Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
        Tχ_upade = get_upade(Tχ,f=f)
        sc = 1 #sqrt(1+a^2+b^2+c^2) # scale to w/J and JS

        for (j,x) in enumerate(x_vec) #J/T
            x0 = x/sc #J1/T

            δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r_max,f=f)[2]

            if x0 == 0.0
                Z[i,j] = 1/Tχ_upade(tanh(f*x0)) * sigma(δ_vec_ext,x0,w*sc,1e-8)
            else
                Z[i,j] = 1/x0 * 1/Tχ_upade(tanh(f*x0)) * sigma(δ_vec_ext,x0,w*sc,1e-8)
            end
        end
    end
    #plt4 = Plots.heatmap(x_vec,a_vec,log10.(Z),xlabel=L"J_x/T",ylabel=L"J_y/J_x",colorbar_title=L"D_S/J_x",aspect_ratio=:auto,c=:balance,shading=:nearest)
    fig = CairoMakie.Figure(size=(325*0.85, 225*0.85))
    ax = CairoMakie.Axis(fig[1, 1],xlabel=L"J_x/T",ylabel=L"J_y/J_x")
    hm = CairoMakie.heatmap!(ax,x_vec,a_vec,Z';colormap=:viridis,colorscale=log10)
    CairoMakie.Colorbar(fig[1, 2],hm;label=L"D_{s,x}/J_x",ticks=[1, 2, 5, 10, 20, 50])
    display(fig)
    #save("Images/SpinDiffusion_$(lattice_type)_heatmap.svg",fig; px_per_unit=6.0)
    #savefig(plt4,"Images/SpinDiffusion_$(lattice_type)_heatmap.png")
end

# comparison D_S/J for different a
function DS_x_0(a_vec::Vector{Float64},str::String;b::Float64=0.0,c::Float64=0.0,r_max=Int64=3,f::Float64=0.7)
    plt5 = Plots.plot([0],[0],label="",xlabel=L"J2/J1",ylabel=L"D_S/J(x=0)",title="Spin Diffusion, $(lattice_type)",size=(1.5*aps_width,aps_width),legend=:bottomright, xlims=(-0.05,4.1)) #xlims=(-0.05,1.1), ylims=(-0.05,0.35)
    
    xs = Float64[]
    ys = Float64[]
    ys_gauß = Float64[]
    for a in a_vec
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
        moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
        c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
        Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
        sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS

        ds_0_gauß = Diffusion_Gauß(Tχ, moments, x/sc)

        if x== 0.0
            ds_0 = 1/Tχ(x) * sigma(get_delta_vec_ext(x,moments,r_max=r_max)[2],x,0.0,0.001)
        else
            ds_0 = 1/x * 1/Tχ(x) * sigma(get_delta_vec_ext(x,moments,r_max=r_max)[2],x,0.0,0.001)
        end

        push!(xs, a)
        push!(ys, ds_0/sc)
        push!(ys_gauß, ds_0_gauß/sc)
    end
    Plots.plot!(plt5, xs, ys, color=:purple, label=false)
    Plots.scatter!(plt5, xs, ys, color=:purple, markersize=6, label="Dyn-HTE")
    Plots.plot!(plt5, xs, ys_gauß, color=:green, label=false)
    Plots.scatter!(plt5, xs, ys_gauß, color=:green, label="Dyn-HTE: Gauß", markershape=:diamond, markersize=6)
    if str == "hypercubic"
        Plots.scatter!(plt5, [0], [sqrt(pi*bc1/(4*d-2-1/(4*bc1)))], label="[Kopietz1993]", markersize=6)
    end
    if lattice_type == "bcc"
        sc2 = sqrt(1+0.6^2)
        x_labrujere = [0.6]
        y_labrujere = [0.155/sc2]
        Plots.scatter!(plt5, x_labrujere, y_labrujere, color=:grey, markershape=:cross, markersize=7, label="Labrujere et al.")
        cap = 0.01
        for (xi, yi, ei) in zip(x_labrujere, y_labrujere, [0.015/sc2])
            plot!([xi, xi], [yi - ei, yi + ei], color=:grey, linewidth=0.8, label=false)
            plot!([xi - cap, xi + cap], [yi + ei, yi + ei], color=:grey, linewidth=0.8, label=false)
            plot!([xi - cap, xi + cap], [yi - ei, yi - ei], color=:grey, linewidth=0.8, label=false)
        end
        
        # add a shaded horizontal band
        y0  = y_labrujere[1]
        half = 0.015/sc2
        #plot!(xs, fill(y0, length(xs)), ribbon=half, color=:grey, alpha=0.2, label="")
        
        #reference data for bcc
        x_M = [0.0]
        y_M = [0.509/2]
        Plots.scatter!(plt5, x_M, y_M, color=:black, markershape=:star5, markersize=9, label="[Morita1972]")
        x_TK = [0.0,0.6]
        y_TK = [0.4916/2,0.5517/(2*sc2)]
        Plots.scatter!(plt5, x_TK, y_TK, color=:orange, markershape=:star5, markersize=9, label="[Tahir-Kheli1970]")
    end
    display(plt5)
    #println(ys)
    #println(ys_gauß)
    #savefig(plt5,"Images/SpinDiffusion_$(lattice_type)_xlim.png")
end

function DS_Bonca(x_vec::Vector{Float64};a::Float64=0.0,b::Float64=0.0,c::Float64=0.0,r_max=Int64=3,f::Float64=0.7)
    plt6 = Plots.plot([0],[0],label="",xlabel=L"T/J",ylabel=L"D_S/J",size=(1.2*aps_width,0.8*aps_width), legend=:topright, ylims=(0.41,0.71))  # ,title="Spin Diffusion, $(lattice_type)_lattice"
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
    moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
    c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
    Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
    sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS

    y_vec = [Diffusion_Gauß(Tχ, moments, x/sc) for x in x_vec]
    Plots.plot!(plt6, 1 ./ x_vec, y_vec, label="Dyn-HTE: Gauß Approx.")

    sig_x = []
    for x in x_vec
        x0 = x/sc
        δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r_max,f=f)[2]

        if x0 == 0.0
            push!(sig_x, 1/Tχ(tanh(f*x0)) * sigma(δ_vec_ext,x0,0.0,0.001))
        else
            push!(sig_x, 1/x0 * 1/Tχ(tanh(f*x0)) * sigma(δ_vec_ext,x0,0.0,0.001))
        end
    end
    
    Plots.plot!(plt6, 1 ./ x_vec, sig_x./sc, label="Dyn-HTE")


    bonca_x = [1.014,1.686,2.340,3.336,5.005,6.672,8.329,10.000]
    bonca_y = [0.642,0.624,0.571,0.522,0.480,0.459,0.447,0.438]
    Plots.scatter!(plt6, bonca_x, bonca_y, yerror=0.05 .* bonca_y, label="[BoncaPRB1995] (FTLM)", color=:grey)

    cap = 0.1
    for (xi, yi, ei) in zip(bonca_x, bonca_y, 0.05 .* bonca_y)
        Plots.plot!([xi, xi], [yi - ei, yi + ei], color=:grey, linewidth=0.8, label=false)
        Plots.plot!([xi - cap, xi + cap], [yi + ei, yi + ei], color=:grey, linewidth=0.8, label=false)
        Plots.plot!([xi - cap, xi + cap], [yi - ei, yi - ei], color=:grey, linewidth=0.8, label=false)
    end

    Plots.plot!(plt6, 1 ./ x_vec, sqrt(pi/20) .* (1 .+ 21/20 .* x_vec), linestyle=:dash, color=:purple, label="[SokolPRL1993]")
    
    # FTLM
    T_vec = [1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 7.5, 8.0, 9.0, 10.0]
    Ds_4x4 = [0.6400959464886398, 0.5856209530002098, 0.5560880450571197, 0.5347289057558636, 0.5067441307966495, 0.48950780344842193, 0.47791165652458495, 0.46959514255775187, 0.4662662597755768, 0.4633507898759228, 0.4584941749009909, 0.4546108039080199]
    Ds_5x5 = [1.2987362361010852, 0.668399892940361, -1.0, 0.5643217912777753, 0.5220882234144044, 0.4987297171627078, 0.4843016994151482, 0.4747443267911697, -1.0, 0.46801563819602, 0.4619956106602265, 0.4581629766630365]
    Plots.scatter!(plt6, T_vec, Ds_4x4, label="[FTLM_4x4]", color=:deepskyblue, saturation=0.5, markershape=:diamond, markersize=4.5)
    Plots.scatter!(plt6, T_vec, Ds_5x5, label="[FTLM_5x5]", color=:blue, markershape=:diamond, markersize=4.5)
    display(plt6)
    #savefig(plt6,"Images/SpinDiffusion_$(lattice_type)_bonca.png")
end

function DS_comp_FTLM(x_vec::Vector{Float64};a::Float64=0.0,b::Float64=0.0,c::Float64=0.0,r_max=Int64=3,f::Float64=0.7)
    plt6 = Plots.plot([0],[0],label="",xlabel=L"J/T",ylabel=L"D_S/J",size=(1.2*aps_width,0.8*aps_width), legend=:topleft, ylims=(0.55,1.85)) # ylims=(0.4,0.7), title="Spin Diffusion, $(lattice_type)_lattice",
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
    moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
    c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
    Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
    Tχ_upade = get_upade(Tχ,f=f)
    sc = 1 # sqrt(1+a^2+b^2+c^2) # scale to w/J and JS

    for r in range(2,4)
        sig_x = []
        for x in x_vec
            x0 = x/sc
            # δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)
            # push!(sig_x, Richardson.extrapolate(w->(1/x0 * 1/Tχ(x0) * sigma_neu(δ_vec,x0,w*sc,extrap_params)),0.5;x0=0.0)[1])
            δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r,f=f)[2]
            if x0 == 0.0
                push!(sig_x, 1/Tχ_upade(tanh(f*x0)) * sigma(δ_vec_ext,x0,0.0,0.001))
            else
                push!(sig_x, 1/x0 * 1/Tχ_upade(tanh(f*x0)) * sigma(δ_vec_ext,x0,0.0,0.001))
            end
        end
        #Plots.plot!(plt6, x_vec, sig_x./sc, label=L"r_\text{max}="*"$r", color=:black, linestyle=linestyle_vec[5-r])
    end

    # FTLM
    T_vec = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    Ds_3x3_T = [2.165734533407567, 1.460535024283905, 1.1653508938529082, 1.1108244321599015, 1.0177603081008768, 1.034914592382753, 0.9843027644121256, 0.9918340170565999, 0.9679016865560089, 0.988302704248027]
    Ds_4x3_T = [1.0012671047847317, 0.8974400637347234, 0.8470310725388679, 0.825580332618162, 0.8296407697720423, 0.7882778255150278, 0.792155904807099, 0.7839655696406969, 0.7938111371315537, 0.7928169447323644]
    Ds_5x3_T = [1.5099417426891455, 0.9753488907505359, 0.8638619346027501, 0.795130203954381, 0.7889863420545117, 0.7909902508317794, 0.776354209484663, 0.7715857504733917, 0.7611591062028444, 0.7562334194174974]
    Ds_4x4_T = [0.6468673995111419, 0.708582693221808, 0.711485985040643, 0.7117027764857609, 0.7116403742851162, 0.7115549602416338, 0.7114822389327528, 0.7114244036398266, 0.711378183303934, 0.7113425481280949]
    Ds_6x3_T = [1.2699266492775947, 0.940006901150526, 0.827479495674566, 0.8153394098593091, 0.7882047309810468, 0.762557474714651, 0.7547022112333623, 0.7535808706212724, 0.7502765742880173, 0.7447624466065175]
    Ds_5x4_T = [1.392346409988413, 0.7883419228983243, 0.7312008867658705, 0.7179269217885581, 0.7095795916114496, 0.7051547849350986, 0.7031728308074237, 0.6986181905756721, 0.7044750953374054, 0.6999801560189247]
    
    b_vec = [0.001, 0.25, 0.5, 0.65, 1.0]
    Ds_4x3_b = [0.7581861366603034, 0.825580332618162, 0.8974400637347234, 0.9295233061386227, 1.0012671047847317]
    Ds_3x4_b = [0.7364308642187003, 0.9224223117910687, 1.1497256863510263, 1.3526167228852561, 2.1115252633835553]
    Ds_4x4_b = [0.7110793461103303 , 0.7117027764857609, 0.708582693221808, 0.7061630462002458, 0.6468673995111419]
    Ds_5x4_b = [0.6968791166910371, 0.7179269217885581, 0.7883419228983243, 0.877889277936424, 1.392346409988413]
    Ds_4x5_b = [0.692015951330529, 0.6970907827530812, 0.7137980462181271, 0.7196066577403811, 0.6877881275574107]
    Ds_5x5_b = [0.6907898838879403, 0.7150205559262233, 0.8011011024238939, 0.8969871215486053, 1.3159447429972093]

    #=
    sat_vec = saturation_palette(6, h=240)
    Plots.plot!(plt6, 1 ./ x_vec, sig_x./sc, label="Dyn-HTE", color=:black)
    # Plots.plot!(plt6,[NaN],[NaN],color=:white,label=" ",linealpha=0,marker=:none)
    Plots.plot!(plt6,[NaN],[NaN],color=:white,label="[FTLM_data]",linealpha=0,marker=:none)
    Plots.scatter!(plt6, T_vec, Ds_3x3_T, label="3x3", color=sat_vec[1], markershape=:diamond, markersize=5)
    Plots.scatter!(plt6, T_vec, Ds_4x3_T, label="4x3", color=sat_vec[2], markershape=:diamond, markersize=5)
    Plots.scatter!(plt6, T_vec, Ds_5x3_T, label="5x3", color=sat_vec[3], markershape=:diamond, markersize=5)
    Plots.scatter!(plt6, T_vec, Ds_4x4_T, label="4x4", color=sat_vec[4], markershape=:diamond, markersize=5)
    Plots.scatter!(plt6, T_vec, Ds_6x3_T, label="6x3", color=sat_vec[5], markershape=:diamond, markersize=5)
    Plots.scatter!(plt6, T_vec, Ds_5x4_T, label="5x4", color=sat_vec[6], markershape=:diamond, markersize=5)
    =#
    #
    sat_vec = saturation_palette(6, h=240)
    #Plots.plot!(plt6, x_vec, sig_x./sc, label="Dyn-HTE", color=:black)
    # Plots.plot!(plt6,[NaN],[NaN],color=:white,label=" ",linealpha=0,marker=:none)
    # Plots.plot!(plt6,[NaN],[NaN],color=:white,label="[FTLM_data]",linealpha=0,marker=:none)
    Plots.scatter!(plt6, b_vec, Ds_4x3_b, label="4x3", color=sat_vec[1], markershape=:diamond, markersize=5)
    # Plots.scatter!(plt6, b_vec, Ds_3x4_b, label="3x4", color=sat_vec[2], markershape=:diamond, markersize=5)
    Plots.scatter!(plt6, b_vec, Ds_4x4_b, label="4x4", color=sat_vec[2], markershape=:diamond, markersize=5)
    Plots.scatter!(plt6, b_vec, Ds_5x4_b, label="5x4", color=sat_vec[3], markershape=:diamond, markersize=5)
    # Plots.scatter!(plt6, b_vec, Ds_4x5_b, label="4x5", color=sat_vec[5], markershape=:diamond, markersize=5)
    Plots.scatter!(plt6, b_vec, Ds_5x5_b, label="5x5", color=sat_vec[6], markershape=:diamond, markersize=5)
    # Plots.hline!([0.9])
    #

    display(plt6)
    #savefig(plt6,"Images/SpinDiffusion_$(lattice_type)_FTLM.svg")

    # inset
    #=
    plt7 = Plots.plot([0],[0],label="",xlabel=L"N",ylabel=L"D_S/J",size=(0.8*aps_width,0.5*aps_width), legend=:topright, xlims=(11,26), ylims=(0.65,1.25), xticks=([12,16,20,25],["12","16","20","25"]))
    Plots.scatter!(plt7,[12,16,20,25],[Ds_4x3_b[1],Ds_4x4_b[1],Ds_5x4_b[1],Ds_5x5_b[1]], color=color_vec, markersize=5, label="J/T=0.0")
    Plots.scatter!(plt7,[12,16,20,25],[Ds_4x3_b[2],Ds_4x4_b[2],Ds_5x4_b[2],Ds_5x5_b[2]], color=color_vec, markershape=:cross, markersize=5, label="J/T=0.25")
    Plots.scatter!(plt7,[12,16,20,25],[Ds_4x3_b[4],Ds_4x4_b[4],Ds_5x4_b[4],Ds_5x5_b[4]], color=color_vec, markershape=:dtriangle, markersize=5, label="J/T=0.65")
    display(plt7)
    # savefig(plt7,"Images/SpinDiffusion_$(lattice_type)_FTLM_inset.png")
    =#
    #=
    plt7 = Plots.plot([0],[0],label="",xlabel=L"N",ylabel=L"D_S/J",size=(0.8*aps_width,0.5*aps_width), legend=:false, xlims=(11,26), ylims=(0.55,1.05), xticks=([12,16,20,25],["12","16","20","25"]))
    Ns = [12, 16, 20, 25]
    Ds_by_N = [Ds_4x3_b, Ds_4x4_b, Ds_5x4_b, Ds_5x5_b]
    idxs = [1, 2, 3, 4]
    labels = ["J/T=0.0", "J/T=0.25", "J/T=0.5", "J/T=0.65"]
    markers = [:circle, :cross, :dtriangle, :star5, :diamond]
    for (idx, label, marker) in zip(idxs, labels, markers)
        yvals = [D[idx] for D in Ds_by_N]
        Plots.scatter!(plt7,Ns,yvals;color=sat_vec,markershape=marker,markersize=5,label=label)
    end
    display(plt7)
    =#
end

# comparison D_S/J for different nmax
function D_nmax(str::String;x::Float64=0.0,f::Float64=0.7,a::Float64=0.0,b::Float64=0.0,c::Float64=0.0)
    if str == "square"
        plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J(x=0)",title="Spin Diffusion, $(lattice_type)",size=(1.5*aps_width,aps_width),legend=:bottomright, xlims=(1.9,6.1),ylims=(0.35,0.5))  # ,title="Spin Diffusion, $(lattice_type)", ylabel=L"D_S/J\;(T\to∞)"
        Plots.scatter!(plt5, [2], [sqrt(pi*bc1/(4*d-2-1/(4*bc1)))], label="[Kopietz1993]", markersize=10)
        data_morita = [0.793,0.842,0.839]
        Plots.scatter!(plt5, [2.0,3.0,4.0], data_morita./2, label="[Morita1972]", markershape=:star, markersize=14)
        Plots.hline!([0.42], color=:red, linestyle=:dash, label=false, linewidth=0.8)
        r_list = [1,2,3,4,5]
    
    elseif str == "simple_cubic"
        plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J(x=0)",title="Spin Diffusion, $(lattice_type)",size=(1.5*aps_width,aps_width),legend=:bottomright, xlims=(1.9,5.1), ylims=(0.25,0.35))
        Plots.scatter!(plt5, [2], [sqrt(pi*bc1/(4*d-2-1/(4*bc1)))], label="[Kopietz1993]", markersize=10)
        data_morita = [0.591,0.620,0.618]
        Plots.scatter!(plt5, [2.0,3.0,4.0], data_morita./2, label="[Morita1972]", markershape=:star, markersize=14)
        Plots.hline!([0.619/2], color=:red, linestyle=:dash, label=false, linewidth=0.8)
        r_list = [1,2,3,4]

    elseif str == "chain"
        plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J(x=0)",title="Spin Diffusion, $(lattice_type)",size=(1.5*aps_width,aps_width),legend=:bottomright, xlims=(1.9,6.1), ylims=(0.75,1.25))
        Plots.scatter!(plt5, [2], [sqrt(pi*bc1/(4*1-2-1/(4*bc1)))], label="[Kopietz1993]", markersize=10)
        data_morita = [1.772,2.105,1.920,2.190]
        Plots.scatter!(plt5, [2.0,3.0,4.0,5.0], data_morita./2, label="[Morita1972]", markershape=:star, markersize=14)
        Plots.hline!([1.05], color=:red, linestyle=:dash, label=false, linewidth=0.8)
        r_list = [1,2,3,4,5]

    elseif str == "bcc"
        plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J(x=0)",title="Spin Diffusion, $(lattice_type)",size=(1.5*aps_width,aps_width),legend=:bottomright, xlims=(1.9,5.1), ylims=(0.15,0.35))
        #Plots.scatter!(plt5, [2], [0.18], label="[Kopietz1993]", markersize=10)
        data_morita = [0.492,0.512,0.507]
        Plots.scatter!(plt5, [2.0,3.0,4.0], data_morita./2, label="[Morita1972]", markershape=:star, markersize=14, color=:green)
        Plots.hline!([0.509/2], color=:red, linestyle=:dash, label=false, linewidth=0.8)
        r_list = [1,2,3,4]

    elseif str == "ladder"
        plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J(x=0)",size=(1.5*aps_width,aps_width),legend=:bottomright, xlims=(1.9,7.1), ylims=(0.75,1.15))
        #Plots.scatter!(plt5, [2], [0.89], label="[Steinigeweg2014]", markersize=8)
        
        # Plot shaded band using ribbon
        xvals = [1.9, 7.1]
        ymid = (0.98 + 0.96) / 2
        yribbon = (0.98 - 0.96) / 2
        Plots.plot!(plt5, xvals, fill(ymid, 2), ribbon=fill(yribbon, 2), fillalpha=0.2, color=:purple, linestyle=:dash, label="[Pollmann2022]")
        
        Plots.hline!([0.95], color=:blue, linestyle=:dash, label="[Reichman2018]", linewidth=0.8)
        r_list = [1,2,3,4,5,6]

    elseif str == "square_XX"
        plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J(x=0)",title="Spin Diffusion, $(lattice_type)",size=(1.5*aps_width,aps_width),legend=:bottomright, xlims=(1.9,7.1))  #ylims=(0.65,0.75),title="Spin Diffusion, $(lattice_type)", ylabel=L"D_S/J\;(T\to∞)"
        r_list = [1,2,3,4]
    end

    xs = Float64[]
    ys = Float64[]
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a)
    moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
    c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
    Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
    Tχ_upade = get_upade(Tχ,f=f)
    sc = 1 #sqrt(1+a^2+b^2+c^2)
    x0 = x/sc

    for r in r_list
        delta_vec_ext = get_delta_vec_ext(x0,moments,r_max=r,f=f)[2]

        if x0 == 0.0
            ds_0 = 1/Tχ_upade(tanh(f*x0)) * sigma(delta_vec_ext,x0,0.0,0.001)
        else
            ds_0 = 1/x0 * 1/Tχ_upade(tanh(f*x0)) * sigma(delta_vec_ext,x0,0.0,0.001)
        end

        push!(xs, r+1)
        push!(ys, ds_0)
    end

    Plots.scatter!(plt5, xs, ys, label="Dyn-HTE", alpha=0.7, color=:black, markersize=6)
    display(plt5)
    println(ys)
    #savefig(plt5,"Images/SpinDiffusion_$(lattice_type)_n.png")
    
    #println(ys)
    #println(data_morita./2)
end

function sus(a::Float64, x_vec::Vector{Float64};f::Float64=0.3,b::Float64=0.0,c::Float64=0.0,m::Int64=6,n::Int64=5)
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
    moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
    c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
    Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:size(c_iipEqualTime_mat)[1]) for n in 0:n_max])  #hte_lattice.lattice.length
    sc = 1 # sqrt(1+a^2+b^2+c^2) # scale to T/J and Jχ

    x_vec_plot = x_vec ./ sc
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
    u_vec = tanh.(f .* x_vec_plot)

    plt = Plots.plot([0],[0],label="",xlabel=L"J/T",ylabel=L"Jχ",size=(0.8*1.5*aps_width,0.8*aps_width), legend=:topleft) 
    #plot!(plt, 1 ./ x_vec_plot, x_vec_plot .* Tχ.(x_vec_plot), label="bare series")
    Plots.plot!(plt, x_vec_plot, x_vec_plot .* get_upade(Tχ,f=f,m=7,n=6).(u_vec), color=color_vec[1], label="Dyn-HTE")
    Plots.plot!(plt, x_vec_plot, x_vec_plot ./ 4, label="J/(4T)", linestyle=:dash, color=color_vec[2])

    ####################
    # Korea
    ####################
    #=
    slope = 1 ./ [454.75, 207.7, 176.18, 124.35, 103.81]
    J = 1/0.47 # in ms
    dV_dx = 2*pi*1.76 # in kHz
    g = 1/0.05  # J / dV_dx
    # println(g)
    for y in slope
        Plots.hline!([g * y], color=:grey, linestyle=:dash, label=false)
    end
    =#

    #display(plt)
    #savefig(plt,"Images/Static_susceptibility.svg")

    # slope
    sites = collect(-10:1:10)
    dV_dx_J = 0.05
    x_vec = [0.0, 0.5, 1.0, 2.0]

    plt2 = Plots.plot([0],[0],label="",xlabel="Position x",ylabel=L"n",size=(0.8*1.5*aps_width,0.8*aps_width), legend=:topleft, ylims=(0.3, 0.7)) 
    colors = cgrad(:coolwarm, length(x_vec), categorical=true)
    for (i,x) in enumerate(x_vec)
        dn_dx = (x * get_upade(Tχ,f=f,m=7,n=6)(tanh(f*x)) * dV_dx_J)
        n = (dn_dx .* sites) .+ 0.5
        Plots.plot!(plt2, sites, n, label="J/T=$x", color=colors[length(x_vec)-i+1])
    end
    display(plt2)
    #savefig(plt2,"Images/dn_dx.svg")
end

####### main ######
x_vec = collect(0.0:0.5:2.0) #J/T
w_vec = collect(0.0:0.025:3.0) #ω/J
a = 1.0
a_vec = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0] #[0.0,0.5,1.1,1.25,1.3,1.4,1.6,2.0,4.0]


##### moments of σ(w) and R(w) are connected but shifted by one -> r_max=3 means we use 0th moment:m_{d,2}, 1st moment:m_{d,4}, 2nd moment m_{d,6}, 3rd moment m_{d,8}

#sigma_w_x([1.2],w_vec,a,r_max=4)
#DS_w_x([0.0],w_vec,a,r_max=6)

#comp_DS_Dyn_Gauß(collect(0.0:0.01:1.0),a,"square",r_max=3)

#sigma_w_a(0.0,exp10.(range(-9, 1.5, length=100)),[0.001,0.005,0.01,0.025,0.05,0.1,0.2],r_max=3)  # [0.001,0.005,0.01,0.025,0.05,0.1,0.2]
#DS_w_a(0.0,w_vec,[0.01,0.1,0.5],r_max=3)
#DS_w_a_moments(0.0,[0.01,0.001,1.0e-6])
#DS_w_a_int(0.0,collect(0.0:1e-3:3.0),[0.1,0.25,0.5],r_max=3)

#DS_x_0([0.0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.5,2.0,2.5,3.0,3.5,4.0],"hypercubic",r_max=3) #[0.0,0.2,0.4,0.5,0.7,0.8,0.9,1.0,1.5,2.0,2.5,3.0,3.5,4.0]
#D_nmax("square_XX",x=1.5,a=1.0)
#sus(1.0, collect(0.0:0.01:1.5))

#DS_Bonca(collect(0.1:0.01:1.0),r_max=3) #only for Heisenberg square lattice
DS_comp_FTLM(collect(0.0:0.01:1.5),r_max=3,a=1.0)

#DS_x_a(collect(0.0:0.01:1.5),0.0,[1.0],r_max=3,f=0.7)
#DS_a_x([1.0,0.5,0.0],0.0,[0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.5,3.0,3.5,4.0,7.0,10.0],r_max=3)
#DS_x_a_heatmap(collect(0.0:0.01:1.0),0.0,[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],r_max=3) # [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]


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
plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J(w=0)",size=(0.8*1.5*aps_width,0.8*aps_width),legend=:bottomright, xlims=(1.9,7.1), ylims=(0.25,1.05)) # title="Spin Diffusion"
Plots.scatter!(plt5, [2,3,4,5,6,7], [0.9594092957640589, 0.995147257146931, 0.9813112720631433, 0.9239983657873552, 0.9402926951404278, 0.9586162740957986], label="ladder", alpha=0.7, color=color_vec, markersize=6)
Plots.scatter!(plt5, [2,3,4,5,6,7], [0.6786602297925832, 0.7039581360308178, 0.7166268745593847, 0.7151146533334913, 0.722074926630639, 0.7146750681412327], label="square lattice", alpha=0.7, color=color_vec, markershape=:cross, markersize=7)
# Plot shaded band using ribbon
xvals = [1.9, 7.1]
ymid = (0.98 + 0.96) / 2
yribbon = (0.98 - 0.96) / 2
Plots.plot!(plt5, xvals, fill(ymid, 2), ribbon=fill(yribbon, 2), fillalpha=0.2, color=color_vec[7], linestyle=:dash, label="[Rakovszky...2022]")
# Plots.hline!([0.95], color=color_vec[8], linestyle=:dash, label="[Reichman2018]", linewidth=0.8)
Plots.hline!([0.72], color=:black, linestyle=:dash, label=false, linewidth=0.8)
display(plt5)
#savefig(plt5,"Images/SpinDiffusion_XX_n.svg")
=#


#=
x_vec = [0.0, 0.5, 1.0, 1.2, 1.5]  # 2.0

b_x0 = [0.6786602297925832, 0.7039581360308178, 0.7166268745593847, 0.7151146533334913]
b_x1 = [0.7321337757622363, 0.7922778626477865, 0.8196388146450287, 0.8192849061332413]
b_x2 = [0.9035213814043781, 1.112771591320647, 1.2146766344036706, 1.2272025093044046]
b_x3 = [1.012310309858948, 1.3473712256480022, 1.5224367956837153, 1.5603613744040197]
# b_x4 = [1.2298474589814432, 1.8947324098460554, 2.291293484056804, 2.460060873806696]
b_x4 = [1.230511333439423, 1.896355725753728, 2.2727280750146295, 2.356906731315863]
b_x5 = [1.7783223178870131, 3.8078406928396764, 4.443267314951987, 7.136832054122836]
b_plot = [b_x0, b_x1, b_x2, b_x3, b_x4]

δ_x0 = [0.6108579656331605, 0.7387803275186121, 0.7015973421780628, 0.7243710163479178]
δ_x1 = [0.6604009648826555, 0.836977225025636, 0.8027656049173182, 0.8299433624997165]
δ_x2 = [0.8202176808189129, 1.1981784534561928, 1.193199865947311, 1.243331866505694]
δ_x3 = [0.9224631866391882, 1.4661362924007164, 1.4995463285702173, 1.582019845374311]
# δ_x4 = [1.1287333190967548, 2.0990674874238238, 2.270377068424058, 2.502742408335409]
δ_x4 = [1.1293690917080146, 2.1009386009659052, 2.249066269434427, 2.3869405719873718]
δ_x5 = [1.6590032227085432, 4.356634723482922, -3.372029407732991, 7.408650750980694]
δ_plot = [δ_x0, δ_x1, δ_x2, δ_x3, δ_x4]

colors = cgrad(:thermal, length(x_vec), categorical=true)
plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J(w=0)",size=(0.8*1.5*aps_width,0.8*aps_width),legend=:bottomright, xlims=(1.8,5.2), ylims=(0.55,2.55))
for (i,x) in enumerate(x_vec)
    lbl_δ = i == 1 ? L"\delta_\text{vec}" : ""
    lbl_b = i == 1 ? L"b_\text{vec}" : ""
    Plots.scatter!(plt5, [2,3,4,5], δ_plot[i], label=lbl_δ, alpha=0.7, color=colors[length(x_vec)-i+1], markersize=6)
    Plots.scatter!(plt5, [2,3,4,5], b_plot[i], label=lbl_b, alpha=0.7, color=colors[length(x_vec)-i+1], markershape=:cross, markersize=7)
    #Plots.scatter!(plt5, [NaN], [NaN],label = L"x = %$x",color = colors[length(x_vec)-i+1],markersize = 6)
end
Plots.hline!(plt5, [pi^0.5/(2*2^0.5)], label="Gaussian approx.", color=:black, linestyle=:dash, linewidth=0.8)
display(plt5)
#savefig(plt5,"Images/SpinDiffusion_square_XX_delta_b_T.svg")
=#


#=
plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel="b_vec_ext",size=(0.8*1.5*aps_width,0.8*aps_width),legend=:bottomright,xlims=(0.9,8.1))
x_vec = [0.0,0.5,1.0,1.2,1.5] # 2.0
b_x0 = [0.3535533905932738, 1.0, 1.7320508075688772, 2.3979157616563596, 3.053270994870246]
b_x1 = [0.3530776607691983, 0.9883658160387082, 1.7580577401778674, 2.4305058651868885, 3.096626711131671]
b_x2 = [0.3514977175005489, 0.9544143287198005, 1.8370172762696895, 2.5175480333193074, 3.206872192030122]
b_x3 = [0.3504906260939432, 0.9351016607958011, 1.884076131286564, 2.561393287487109, 3.263052438486597]
b_x4 = [0.34853407021911703, 0.9005898072216952, 1.972325659868817, 2.628812406223184, 3.359442804653361]
b_x5 = [0.3440601606326963, 0.8300561033617464, 2.171202557440214, 2.7194255396667972, 3.570280184201495]
b = [b_x0, b_x1, b_x2, b_x3, b_x4]
colors = cgrad(:thermal, length(x_vec), categorical=true)
for (i,x) in enumerate(x_vec)
    Plots.scatter!(plt5, [1,2,3,4,5], b[i], label="x=$x", color=colors[length(x_vec)-i+1])
    Plots.plot!(plt5, 5:8, extrapolate_δvec(b[i],length(b[i])-2,length(b[i])-1,2000,false)[5:8], color=colors[length(x_vec)-i+1], label=false)
end

display(plt5)
#savefig(plt5,"Images/SpinDiffusion_b_T.svg")
=#


#c_iipDyn_mat_subst = load_ciipDyn_mat_subst(0.0)
#moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
#get_delta_vec_ext_neu(0.0,moments,r_max=5,r_min=3)


#=
c_iipDyn_mat_subst = load_ciipDyn_mat_subst(0.0)
moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
moments_morita = [moments[i] * 2^(2*i) for i in eachindex(moments)]
#println([moments[i](0.0) for i in eachindex(moments)])
println(gauss_moments_ratio(moments,0.0))
Plots.scatter([2,3,4,5,6],gauss_moments_ratio(moments,0.0),xlabel="r",ylabel=L"m_{d,2r}/m_{d,2r}^{\text{Gauss}}", markersize=6,label="",size=(1.5*aps_width*0.5,aps_width*0.5)) # title="moments ratio, $(lattice_type) lattice",
display(current())
Plots.savefig("Images/Moments_Ratio_$(lattice_type).png")
=#







##############################################################################
##############################################################################
##############################################################################
if false
x_vec = 0:0.01:2.0
c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a)
m_vec = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
m_vec_times_x = [m_vec[i]*Polynomial([0,1]) for i=1:length(m_vec)]

##EXTRAPLOATION OF MOMENTS
#basic pade 
m_vec_extrapolated_pade = []
for m_idx=1:length(m_vec) #-2
    push!(m_vec_extrapolated_pade, extrapolate_series(m_vec[m_idx],"pade",(7-m_idx,7-m_idx)))
end

#pade in u = tanh(f*x) (2 different versions)
f = 0.7   #define the f value (f=0.48 is very fine tuned to give good results)
m_vec_times_x_extrapolated_u_pade1 = []
m_vec_times_x_extrapolated_u_pade2 = []
for m_idx=1:length(m_vec) #-2
    push!(m_vec_times_x_extrapolated_u_pade1, extrapolate_series(m_vec_times_x[m_idx],"u_pade",(8-m_idx,7-m_idx,f)))
    push!(m_vec_times_x_extrapolated_u_pade2, extrapolate_series(m_vec_times_x[m_idx],"u_pade",(7-m_idx,8-m_idx,f)))
end

#plot the moments 
plt_m = Plots.plot([0],[0],label="",xlabel=L"x",ylabel=L"x \cdot m_{\mathbf{k},r}(x)/m_{\mathbf{k},r}(0)",title="moments",legend=:topleft,xlim=(-0.2,2.2),size=(1.5*aps_width,aps_width))
Plots.plot!(plt_m,[0],[0],label="x bare",color = "grey",linestyle = linestyle_vec[1],linewidth=0.4)
Plots.plot!(plt_m,[0],[0],label="u Padé [7-r,6-r]",color = "grey",linestyle = linestyle_vec[2],alpha =0.5)
Plots.plot!(plt_m,[0],[0],label="u Padé [6-r,7-r]",color = "grey",linestyle = linestyle_vec[3])
#Plots.plot!(plt_m,[0],[0],label="Padé [6-r,6-r]",color = "grey",linestyle = linestyle_vec[4])
for i=1:4
    Plots.plot!(plt_m,[0],[0],label="r="*string(i-1),color = color_vec[i])
    Plots.plot!(plt_m,x_vec[1:180],m_vec_times_x[i].(x_vec[1:180])./(m_vec[i](0)),label = "",alpha= 0.7,color = color_vec[i],linestyle = linestyle_vec[1],linewidth=0.5)
    Plots.plot!(plt_m,x_vec,m_vec_times_x_extrapolated_u_pade1[i].(tanh.(f.*x_vec))./(m_vec[i](0)),label = "",alpha= 0.5,color = color_vec[i],linestyle = linestyle_vec[2],linewidth=0.5)
    Plots.plot!(plt_m,x_vec,m_vec_times_x_extrapolated_u_pade2[i].(tanh.(f.*x_vec))./(m_vec[i](0)),label = "",alpha= 1,color = color_vec[i],linestyle = linestyle_vec[3],linewidth=0.5)
    #Plots.plot!(plt_m,x_vec,x_vec.*m_vec_extrapolated_pade[i].(x_vec)./(m_vec_extrapolated_pade[i](0.0001)),label = "",alpha= 1,color = color_vec[i],linestyle = linestyle_vec[4],linewidth=0.5)
end

#SHOW PLOT
##############
display(plt_m)

end

###### check deltas
if false

#now calculate the δs of the continued fraction expansion 
#(we will use the u_pade extrapolation)
x_vec = [0.5,1.0,2.0,4.0]  #define the x values at which we calculate the deltas

deltas_at_x =[]
for x in x_vec
    δ_vec,r_vec = fromMomentsToδ([m(tanh(f*x))/x for m in m_vec_times_x_extrapolated_u_pade1])
    push!(deltas_at_x,δ_vec)
end

#plots the δs 
plt_δ = Plots.plot([0],[0],title="δs using Padé with substitution on moments",label="",xlabel=L"r",ylabel=L"\delta_r",legend=:topleft,size=(1.5*aps_width,aps_width))
for (i,x) in enumerate(x_vec)
    Plots.scatter!(plt_δ,0:3,deltas_at_x[i][1:4],label = "x="*string(x),color = thermalCol4_vec[i],markersize=7)
end

#SHOW PLOT
##############
display(plt_δ)

#one can now extrapolate the deltas to infinity via
deltas_at_x_ext = []
for x_idx = 1:length(x_vec)
    δ_vec_ext = extrapolate_δvec(deltas_at_x[x_idx],length(δ_vec)-1,length(δ_vec)-1,2000,false)
    push!(deltas_at_x_ext,δ_vec_ext)
end

#add the extrapolated deltas to the plot 
for (i,x) in enumerate(x_vec)
    Plots.scatter!(plt_δ,4:10,deltas_at_x_ext[i][5:11],label = "",title="δs extrapolated",color = thermalCol4_vec[i],markersize=4,alpha=0.5)
end

#SHOW PLOT
##############
display(plt_δ)

#we can use the deltas to calculate the spin structure factor (still at fixed k)
w_vec = 0.01:0.01:3  #define w
x_idx = 3          #choose temperature index (relativ to x_vec = [0.5,1.0,2.0,4.0]) 
η = 0.02          #the imaginary part after analytic continuation (η ->0)
#now calculate the DSF 
DSF = [sigma(deltas_at_x_ext[x_idx],x_vec[x_idx],w,η) for w in w_vec]

plt_dsf = Plots.plot(w_vec,DSF,label="",title="Spin Conductivity at x="*string(x_vec[x_idx]),xlabel=L"w",ylabel=L"\sigma_S^{xx}(\mathbf{k}=0,w)",size = (1.5*aps_width,aps_width))

#SHOW PLOT
##############
display(plt_dsf)
end


######## old functions #######
#=
# comparison D_S with Gauß and literature
function comp_DS_Dyn_Gauß_old(x_vec::Vector{Float64},w::Float64,a::Float64,str::String;b::Float64=0.0,c::Float64=0.0,r_max=Int64=3,f::Float64=0.7)
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
    moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
    c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
    Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
    sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS

    y_vec = [Diffusion_Gauß(Tχ, moments, x/sc) for x in x_vec]
    plt3 = Plots.plot([0],[0],label="",xlabel=L"J/T",ylabel=L"D_S/J",title="Spin Diffusion, $(lattice_type), J2/J1=$(a)",size=(1.5*aps_width,aps_width), legend=:bottomleft)
    Plots.plot!(plt3, x_vec, y_vec, label="Dyn-HTE: Gauß Approximation")

    #sig_x = [1/(i/sc) * 1/Tχ(i/sc) * sigma(get_delta_vec_ext(i/sc,moments,r_max=r_max,f=f),i/sc,w,0.02) for i in x_vec]
    
    sig_x = []
    for x in x_vec
        x0 = x/sc
        δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)
        push!(sig_x, 1/x0 * 1/Tχ(x0) * sigma_neu(δ_vec,x0,w*sc,extrap_params))
    end
    
    Plots.plot!(plt3, x_vec, sig_x./sc, label="Dyn-HTE: w=$(w)")
    println(sig_x[1]-y_vec[1])

    #### hypercubic lattice
    if str == "hypercubic" && a==0.0
        Plots.scatter!(plt3, [0], [sqrt(pi*bc1/(4*d-2-1/(4*bc1)))], label="[Kopietz1993]", markersize=6)
    elseif str == "square" && a==0.0
    #### square
        Plots.scatter!(plt3, [0], [sqrt(pi*bc1/(4*d-2-1/(4*bc1)))], label="[Kopietz1993]", markersize=6)
        Plots.plot!(plt3,x_vec, sqrt(pi/20) .* (1 .+ 21/20 .* x_vec), linestyle=:dash, label="[SokolPRL1993]")
    else
        println("No literature values")
    end
    display(plt3)
    #savefig(plt3,"Images/SpinDiffusion_$(lattice_type)_lit_a_$(a).png")
end
#comp_DS_Dyn_Gauß_old(collect(0.001:0.01:1.0),wlim,a,"hypercubic",r_max=3)

# comparison D_S(x) for different a
function DS_x_0_old(x::Float64,w::Float64,a_vec::Vector{Float64},str::String;b::Float64=0.0,c::Float64=0.0,r_max=Int64=3,f::Float64=0.7)
    plt5 = Plots.plot([0],[0],label="",xlabel=L"J2/J1",ylabel=L"D_S/J(x=0)",title="Spin Diffusion, $(lattice_type), w=$(w)",size=(1.5*aps_width,aps_width),legend=:bottomright)
    xs = Float64[]
    ys = Float64[]
    for a in a_vec
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
        moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
        c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
        Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
        sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS

        #sig_x = 1/(x/sc) * 1/Tχ(x/sc) * sigma(get_delta_vec_ext(x/sc,moments,r_max=r_max,f=f),x/sc,w,0.02)

        x0 = x/sc
        δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)
        sig_x = 1/x0 * 1/Tχ(x0) * sigma_neu(δ_vec,x0,w*sc,extrap_params)

        push!(xs, a)
        push!(ys, sig_x/sc)
    end
    Plots.plot!(plt5, xs, ys, label=false)
    Plots.scatter!(plt5, xs, ys, label="Dyn-HTE")
    if str == "hypercubic"
        Plots.scatter!(plt5, [0], [sqrt(pi*bc1/(4*d-2-1/(4*bc1)))], label="[Kopietz1993]", markersize=6)
    end
    display(plt5)
    #savefig(plt5,"Images/SpinDiffusion_$(lattice_type)_xlim.png")
end
#DS_x_0_old(xlim,wlim,[0.0,0.2,0.4,0.5,0.7,0.8,0.9,1.0,1.5,2.0,2.5,3.0,3.5,4.0],"hypercubic",r_max=3)
=#