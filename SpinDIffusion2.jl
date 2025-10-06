using JLD2, DelimitedFiles, SimpleWeightedGraphs, Plots, Symbolics, Richardson

include("plotConventions.jl")
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 

@variables x1 x2 Δ

### load graph evaluations
spin_length = 1/2
bc1 = spin_length*(spin_length+1)/3
n_max = 10

### prepare lattice
lattice_type = "bcc" #"Shastry-Sutherland"
d = 3
ex = zeros(d)  # creates a d-dimensional zero vector
ex[1] = 1      # set first component to 1

# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = true
j3 = false
j4 = false

L = 10

#hte_lattice = getLattice(L,lattice_type,j1,j2,j3,j4);

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

function get_moments_Gklim(c_iipDyn_mat, hte_lattice, ex)
    lattice = hte_lattice.lattice
    center_sites = hte_lattice.basis_positions
    z = zeros(size(c_iipDyn_mat[1]))
    for b in 1:length(center_sites)
        for i in 1:length(lattice)
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
    res = 1/pi * real(contFrac(1im * w + η,δ_vec))
    if x==0.0 || w==0.0
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

    res = 1/pi * real(contFracwithTerminator(1im * w ,δ_vec,a,b))
    if x==0.0 || w==0.0
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

function get_delta_vec_ext(x0::Float64, moments;f::Float64=0.7,r_max::Int64=3)
    n_max_d = n_max#-2
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max_d+1,f)
    poly_x = Polynomial([0,1],:x)
    u0 = tanh.(f .* x0)
        m0 = Float64[]
        for r in 1:r_max
            xm_norm_r = coeffs(poly_x * (moments[1+r]/moments[1+r](0)))
            p_u = Polynomial(ufromx_mat[1:n_max_d+2-2*r,1:n_max_d+2-2*r]*xm_norm_r)
            push!(m0,moments[1+r](0)/x0 * get_pade(p_u,7-r,6-r)(u0))
        end
    δ_vec,r_vec = fromMomentsToδ(m0)
    δ_vec_ext = extrapolate_δvec(δ_vec,length(δ_vec)-1,length(δ_vec)-1,2000,true)
    return δ_vec_ext
end

function get_delta_vec_ext_neu(x0::Float64,moments;f::Float64=0.7,r_max::Int64=3,r_min::Int64=3)
    intercept0 = true
    n_max_d = n_max - 2
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max_d+1,f)
    poly_x = Polynomial([0,1],:x)
    u0 = tanh.(f .* x0)
        m0 = Float64[]
        for r in 0:r_max
            xm_norm_r = coeffs(poly_x * (moments[1+r]/moments[1+r](0)))
            p_u = Polynomial(ufromx_mat[1:n_max_d+2-2*r,1:n_max_d+2-2*r]*xm_norm_r)
            push!(m0,moments[1+r](0)/x0 * get_pade(p_u,6-r,5-r)(u0))
        end
    δ_vec,r_vec = fromMomentsToδ(m0)
    
    # find last index where δ_vec is non-negative
    idx = findfirst(<(0), δ_vec)
    lastidx = isnothing(idx) ? length(δ_vec) : idx - 1
    r_max_eff = min(lastidx - 1, r_max)
    r_min_eff = min(r_min, r_max_eff)

    if r_max_eff < 1
        println("WARNING: negative δ1, putting δ0 = 0")
        extrap_params = [1.0,0.0]
        δ_vec = [0.0,1.0]
        r_max_eff = 0
    else
        extrap_params = get_extrapolation_params(δ_vec[1:r_max_eff+1],r_min_eff,r_max_eff,intercept0)
    end
    return δ_vec, extrap_params
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

        #δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r_max,f=f)
        #sigma_ω = [sigma(δ_vec_ext,x0,w*sc,0.01) for w in w_vec]

        δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)
        sigma_ω = [sigma_neu(δ_vec,x0,w*sc,extrap_params) for w in w_vec]

        Plots.plot!(plt, w_vec, sigma_ω, label="x=$(x)")
    end
    display(plt)
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
    
    plt = Plots.plot([0],[0],label="",xlabel=L"w",ylabel=L"D_S/J(w)",title="$(lattice_type), J2/J1=$(a)",size=(1.5*aps_width,aps_width), legend=:topright)

    for x in x_vec # J/T
        x0 = x/sc  # J1/T

        #δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r_max,f=f)
        #sigma_ω = [sigma(δ_vec_ext,x0,w*sc,0.01) for w in w_vec]

        δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)

        #DS_w = [1/x0 * Tχ(x0) * sigma_neu(δ_vec,x0,w*sc,extrap_params) for w in w_vec]
        DS_w = [1/x0 * Tχ_upade(tanh(f*x0)) * sigma_neu(δ_vec,x0,w*sc,extrap_params) for w in w_vec]

        Plots.plot!(plt, w_vec, DS_w, label="x=$(x)")
    end
    display(plt)
    #savefig(plt,"Images/DS_$(lattice_type)_a_$(a).png")
end

function DS_w_x_log(x_vec::Vector{Float64},w_vec::Vector{Float64},a::Float64;r_max::Int64=3,f::Float64=0.7,b::Float64=0.0,c::Float64=0.0)
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
    moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
    c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
    Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
    Tχ_upade = get_upade(Tχ,f=f)

    sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS
    
    plt = Plots.plot([0],[0],label="",xlabel=L"w^2",ylabel=L"D_S/J(w)",title="$(lattice_type), J2/J1=$(a)",size=(1.5*aps_width,aps_width), legend=:topright, ylim=(0.0001,1.0), yscale=:log)

    for x in x_vec # J/T
        x0 = x/sc  # J1/T

        #δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r_max,f=f)
        #sigma_ω = [sigma(δ_vec_ext,x0,w*sc,0.01) for w in w_vec]

        δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)

        #DS_w = [1/x0 * Tχ(x0) * sigma_neu(δ_vec,x0,w*sc,extrap_params) for w in w_vec]
        DS_w = [1/x0 * Tχ_upade(tanh(f*x0)) * sigma_neu(δ_vec,x0,w*sc,extrap_params) for w in w_vec]

        Plots.plot!(plt, w_vec.^2, DS_w, label="x=$(x)")
    end
    display(plt)
    #savefig(plt,"Images/DS_$(lattice_type)_a_$(a)_log.png")
end

# comparison comparison σ(w) for different values of a at fixed x
function sigma_w_a(x::Float64, w_vec::Vector{Float64}, a_vec::Vector{Float64};b::Float64=0.0,c::Float64=0.0,r_max::Int64=3,f::Float64=0.7)
    plt2 = Plots.plot([0],[0],label="",xlabel=L"w",ylabel=L"\sigma(w)",title="Spin conductivity, $(lattice_type), x=$(x)",size=(1.5*aps_width,aps_width), legend=:topright)
    for a in a_vec
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
        moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
        sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS
        x0 = x/sc # J1/T
        
        #δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r_max,f=f)
        #sigma_ω = [sigma(δ_vec_ext,x0,w*sc,0.01) for w in w_vec]

        δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)
        sigma_ω = [sigma_neu(δ_vec,x0,w*sc,extrap_params) for w in w_vec]

        Plots.plot!(plt2, w_vec, sigma_ω, label="a=$(a)")
        #println(sigma_ω[1]/(Tχ(x0)*x0))
    end
    display(plt2)
    #savefig(plt2,"Images/SpinConductivity_$(lattice_type)_x_$(x).png")
end

# comparison comparison D_S/J(w) for different values of a at fixed x
function DS_w_a(x::Float64, w_vec::Vector{Float64}, a_vec::Vector{Float64};b::Float64=0.0,c::Float64=0.0,r_max::Int64=3,f::Float64=0.7)
    plt2 = Plots.plot([0],[0],label="",xlabel=L"w",ylabel=L"D_S/J(w)",title="$(lattice_type), x=$(x)",size=(1.5*aps_width,aps_width), legend=:topright)
    for a in a_vec
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
        moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
        c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
        Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
        Tχ_upade = get_upade(Tχ,f=f)
        sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS
        x0 = x/sc # J1/T
        
        #δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r_max,f=f)
        #sigma_ω = [sigma(δ_vec_ext,x0,w*sc,0.01) for w in w_vec]

        δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)
        
        #sigma_ω = [1/x0 * Tχ(x0) * sigma_neu(δ_vec,x0,w*sc,extrap_params) for w in w_vec]
        sigma_ω = [1/x0 * Tχ_upade(tanh(f*x0)) * sigma_neu(δ_vec,x0,w*sc,extrap_params) for w in w_vec]

        Plots.plot!(plt2, w_vec, sigma_ω, label="a=$(a)")
        #println(sigma_ω[1]/(Tχ(x0)*x0))
    end
    display(plt2)
    #savefig(plt2,"Images/DS_$(lattice_type)_x_$(x).png")
end

# comparison comparison D_S/J(w) for different values of a at fixed x
function DS_w_a_log(x::Float64, w_vec::Vector{Float64}, a_vec::Vector{Float64};b::Float64=0.0,c::Float64=0.0,r_max::Int64=3,f::Float64=0.7)
    plt2 = Plots.plot([0],[0],label="",xlabel=L"w^2",ylabel=L"D_S/J(w)",title="$(lattice_type), x=$(x)",size=(1.5*aps_width,aps_width), legend=:topright, ylim=(0.0001,1.0), yscale=:log)
    for a in a_vec
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
        moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
        c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
        Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
        Tχ_upade = get_upade(Tχ,f=f)
        
        sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS
        x0 = x/sc # J1/T
        
        #δ_vec_ext = get_delta_vec_ext(x0,moments,r_max=r_max,f=f)
        #sigma_ω = [sigma(δ_vec_ext,x0,w*sc,0.01) for w in w_vec]

        δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)

        #sigma_ω = [1/x0 * Tχ(x0) * sigma_neu(δ_vec,x0,w*sc,extrap_params) for w in w_vec]
        sigma_ω = [1/x0 * Tχ_upade(tanh(f*x0)) * sigma_neu(δ_vec,x0,w*sc,extrap_params) for w in w_vec]

        Plots.plot!(plt2, w_vec .^ 2, sigma_ω, label="a=$(a)")
    end
    display(plt2)
    #savefig(plt2,"Images/DS_$(lattice_type)_x_$(x)_log.png")
end

# comparison D_S with Gauß and literature
function comp_DS_Dyn_Gauß(x_vec::Vector{Float64},a::Float64,str::String;b::Float64=0.0,c::Float64=0.0,r_max=Int64=3,f::Float64=0.7)
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
    moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
    c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
    Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
    Tχ_upade = get_upade(Tχ,f=f)
    sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS

    y_vec = [Diffusion_Gauß(Tχ, moments, x/sc) for x in x_vec]
    plt3 = Plots.plot([0],[0],label="",xlabel=L"J/T",ylabel=L"D_S/J",title="Spin Diffusion, $(lattice_type), J2/J1=$(a)",size=(1.5*aps_width,aps_width), legend=:bottomleft)
    Plots.plot!(plt3, x_vec, y_vec, label="Dyn-HTE: Gauß Approximation")

    #sig_x = [1/(i/sc) * 1/Tχ(i/sc) * sigma(get_delta_vec_ext(i/sc,moments,r_max=r_max,f=f),i/sc,w,0.02) for i in x_vec]
    
    sig_x = []
    for x in x_vec #J/T
        x0 = x/sc # J1/T
        δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)
        
        push!(sig_x, 1/x0 * 1/Tχ(x0) * Richardson.extrapolate(w->sigma_neu(δ_vec,x0,w*sc,extrap_params),0.5;x0=0.0)[1])
        #push!(sig_x, 1/x0 * 1/Tχ_upade(tanh(f*x0)) * Richardson.extrapolate(w->sigma_neu(δ_vec,x0,w*sc,extrap_params),0.5;x0=0.0)[1])
        
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
    plt4 = Plots.plot([0],[0],label="",xlabel=L"J/T",ylabel=L"D_S/J",title="Spin Diffusion, $(lattice_type), w=$(w)",size=(1.5*aps_width,aps_width), legend=:topleft)
    for a in a_vec
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
        moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
        c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
        Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
        Tχ_upade = get_upade(Tχ,f=f)
        sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS

        #sig_x = [1/(i/sc) * 1/Tχ(i/sc) * sigma(get_delta_vec_ext(i/sc,moments,r_max=r_max,f=f),i/sc,w,0.02) for i in x_vec]
        
        sig_x = []
        for x in x_vec #J/T
            x0 = x/sc #J1/T
            δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)
            
            #push!(sig_x, 1/x0 * 1/Tχ(x0) * sigma_neu(δ_vec,x0,w*sc,extrap_params))
            push!(sig_x, 1/x0 * 1/Tχ_upade(tanh(f*x0)) * sigma_neu(δ_vec,x0,w*sc,extrap_params))
        end
        Plots.plot!(plt4, x_vec, sig_x./sc, label="J2/J1=$(a)")
    end
    display(plt4)
    #savefig(plt4,"Images/SpinDiffusion_$(lattice_type).png")
end

# comparison D_S/J for different a
function DS_x_0(a_vec::Vector{Float64},str::String;b::Float64=0.0,c::Float64=0.0,r_max=Int64=3,f::Float64=0.7)
    plt5 = Plots.plot([0],[0],label="",xlabel=L"J2/J1",ylabel=L"D_S/J(x=0)",title="Spin Diffusion, $(lattice_type)",size=(1.5*aps_width,aps_width),legend=:bottomright, xlims=(-0.05,1.1), ylims=(-0.05,0.35)) #xlims=(-0.05,1.1), ylims=(-0.05,0.35)
    
    xs = Float64[]
    ys = Float64[]
    ys_gauß = Float64[]
    for a in a_vec
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
        moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
        c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
        Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
        sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS

        ds_0_gauß = Richardson.extrapolate(x->Diffusion_Gauß(Tχ, moments, x/sc),0.5;x0=0.0)[1]
        ds_0 = Richardson.extrapolate(x->(1/x * 1/Tχ(x) * Richardson.extrapolate(w->sigma_neu(get_delta_vec_ext_neu(x,moments,r_max=r_max,f=f)[1],x,w*sc,get_delta_vec_ext_neu(x,moments,r_max=r_max,f=f)[2]),0.5;x0=0.0)[1]),0.1;x0=0.0)[1]

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
    plt6 = Plots.plot([0],[0],label="",xlabel=L"T/J",ylabel=L"D_S/J",title="Spin Diffusion, $(lattice_type)_lattice",size=(1.5*aps_width,aps_width), legend=:bottomleft, ylims=(0.4,0.7))
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
    moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
    c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
    Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
    sc = sqrt(1+a^2+b^2+c^2) # scale to w/J and JS

    y_vec = [Diffusion_Gauß(Tχ, moments, x/sc) for x in x_vec]
    Plots.plot!(plt6, 1 ./ x_vec, y_vec, label="Dyn-HTE: Gauß Approx.")

    #sig_x = [1/(i/sc) * 1/Tχ(i/sc) * sigma(get_delta_vec_ext(i/sc,moments,r_max=r_max,f=f),i/sc,w,0.02) for i in x_vec]
    
    sig_x = []
    for x in x_vec
        x0 = x/sc
        δ_vec,extrap_params = get_delta_vec_ext_neu(x0,moments,r_max=r_max,f=f)
        push!(sig_x, Richardson.extrapolate(w->(1/x0 * 1/Tχ(x0) * sigma_neu(δ_vec,x0,w*sc,extrap_params)),0.5;x0=0.0)[1])
    end
    
    Plots.plot!(plt6, 1 ./ x_vec, sig_x./sc, label="Dyn-HTE")

    bonca_x = [1.014,1.686,2.340,3.336,5.005,6.672,8.329,10.000]
    bonca_y = [0.642,0.624,0.571,0.522,0.480,0.459,0.447,0.438]
    Plots.scatter!(plt6, bonca_x, bonca_y, yerror=0.05 .* bonca_y, label="[BoncaPRB1995] (FTLM)", color=:grey)

    cap = 0.1
    for (xi, yi, ei) in zip(bonca_x, bonca_y, 0.05 .* bonca_y)
        plot!([xi, xi], [yi - ei, yi + ei], color=:grey, linewidth=0.8, label=false)
        plot!([xi - cap, xi + cap], [yi + ei, yi + ei], color=:grey, linewidth=0.8, label=false)
        plot!([xi - cap, xi + cap], [yi - ei, yi - ei], color=:grey, linewidth=0.8, label=false)
    end

    Plots.plot!(plt6, 1 ./ x_vec, sqrt(pi/20) .* (1 .+ 21/20 .* x_vec), linestyle=:dash, color=:purple, label="[SokolPRL1993]")
    
    display(plt6)
    #savefig(plt6,"Images/SpinDiffusion_$(lattice_type)_bonca.png")
end

# comparison D_S/J for different nmax
function D_nmax(str::String;f::Float64=0.7)
    if str == "square"
        plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J\;(T\to∞)",size=(1.5*aps_width,aps_width),legend=:bottomright, xlims=(1.9,5.1),ylims=(0.39,0.45))  # ,title="Spin Diffusion, $(lattice_type)"
        Plots.scatter!(plt5, [2], [sqrt(pi*bc1/(4*d-2-1/(4*bc1)))], label="[Kopietz1993]", markersize=10)
        data_morita = [0.793,0.842,0.839]
        Plots.scatter!(plt5, [2.0,3.0,4.0], data_morita./2, label="[Morita1972]", markershape=:star, markersize=14)
        Plots.hline!([0.42], color=:red, linestyle=:dash, label=false, linewidth=0.8)
        r_list = [1,2,3,4]
    
    elseif str == "sc"
        plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J(x=0)",title="Spin Diffusion, $(lattice_type)",size=(1.5*aps_width,aps_width),legend=:bottomright, xlims=(1.9,5.1), ylims=(0.25,0.35))
        Plots.scatter!(plt5, [2], [sqrt(pi*bc1/(4*d-2-1/(4*bc1)))], label="[Kopietz1993]")
        data_morita = [0.591,0.620,0.618]
        Plots.scatter!(plt5, [2.0,3.0,4.0], data_morita./2, label="[Morita1972]", markershape=:star, markersize=8)
        Plots.hline!([0.619/2], color=:red, linestyle=:dash, label=false, linewidth=0.6)
        r_list = [1,2,3]

    elseif str == "chain"
        plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J(x=0)",title="Spin Diffusion, $(lattice_type)",size=(1.5*aps_width,aps_width),legend=:bottomright, xlims=(1.9,5.1), ylims=(0.75,1.25))
        Plots.scatter!(plt5, [2], [sqrt(pi*bc1/(4*1-2-1/(4*bc1)))], label="[Kopietz1993]")
        data_morita = [1.772,2.105,1.920,2.190]
        Plots.scatter!(plt5, [2.0,3.0,4.0,5.0], data_morita./2, label="[Morita1972]", markershape=:star, markersize=8)
        Plots.hline!([1.05], color=:red, linestyle=:dash, label=false, linewidth=0.6)
        r_list = [1,2,3,4]

    elseif str == "bcc"
        plt5 = Plots.plot([0],[0],label="",xlabel=L"r_{max}",ylabel=L"D_S/J(x=0)",title="Spin Diffusion, $(lattice_type)",size=(1.5*aps_width,aps_width),legend=:bottomright, xlims=(1.9,5.1), ylims=(0.15,0.35))
        #Plots.scatter!(plt5, [2], [0.18], label="[Kopietz1993]")
        data_morita = [0.492,0.512,0.507]
        Plots.scatter!(plt5, [2.0,3.0,4.0], data_morita./2, label="[Morita1972]", markershape=:star, markersize=8, color=:green)
        Plots.hline!([0.509/2], color=:red, linestyle=:dash, label=false, linewidth=0.6)
        r_list = [1,2,3]
    end

    xs = Float64[]
    ys = Float64[]
    for r in r_list
        c_iipDyn_mat_subst = load_ciipDyn_mat_subst(0.0)
        moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
        c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
        Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])

        ds_0 = Richardson.extrapolate(x->(1/x * 1/Tχ(x) * Richardson.extrapolate(w->sigma_neu(get_delta_vec_ext_neu(x,moments,r_max=r,f=f)[1],x,w,get_delta_vec_ext_neu(x,moments,r_max=r,f=f)[2]),0.5;x0=0.0)[1]),0.1;x0=0.0)[1]
        push!(xs, r+1)
        push!(ys, ds_0)
    end
    Plots.scatter!(plt5, xs, ys, label="Dyn-HTE", alpha=0.7, color=:black, markersize=6)
    display(plt5)
    #savefig(plt5,"Images/SpinDiffusion_$(lattice_type)_n.png")
    
    #println(ys)
    #println(data_morita./2)
end

function sus(a::Float64, x_vec::Vector{Float64};f::Float64=0.3,b::Float64=0.0,c::Float64=0.0,m::Int64=6,n::Int64=5)
    c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a,b=b,c=c)
    moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
    c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
    Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
    sc = sqrt(1+a^2+b^2+c^2) # scale to T/J and Jχ

    x_vec_plot = x_vec ./ sc
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
    u_vec = tanh.(f .* x_vec_plot)

    plt = Plots.plot([0],[0],label="",xlabel=L"T/J",ylabel=L"Jχ",title="$(lattice_type), a=$(a)",size=(1.5*aps_width,aps_width), legend=:topright) 
    #plot!(plt, 1 ./ x_vec_plot, x_vec_plot .* Tχ.(x_vec_plot), label="bare series")
    plot!(plt, 1 ./ x_vec_plot, x_vec_plot .* get_upade(Tχ,f=f,m=m,n=n).(u_vec), label="u-pade[$(m),$(n)]")

    #scatter!(plt,[0.5,1.0,1.5,2.0],get_upade(Tχ,f=f,m=5,n=4).(tanh.(f.*[0.5,1.0,1.5,2.0])))
    display(plt)
    #println(Tχ)
end

####### main ######
x_vec = collect(0.01:0.5:1.5) #J/T
w_vec = collect(0.025:0.025:5.0) #ω/J
a = 0.0
a_vec = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0] #[0.0,0.5,1.1,1.25,1.3,1.4,1.6,2.0,4.0]


# moments of σ(w) and R(w) are connected but shifted by one -> r_max=3 means we use 0th moment:m_{d,2}, 1st moment:m_{d,4}, 2nd moment m_{d,6}, 3rd moment m_{d,8}

#sigma_w_x(x_vec,w_vec,a,r_max=3)
#DS_w_x(x_vec,w_vec,a,r_max=3)
#DS_w_x_log(x_vec,w_vec,a,r_max=3)

#comp_DS_Dyn_Gauß(collect(0.001:0.01:1.0),a,"none",r_max=3)

#sigma_w_a(0.01,w_vec,a_vec,r_max=3)
#DS_w_a(1.0,w_vec,a_vec,r_max=3)
#DS_w_a_log(0.01,w_vec,a_vec,r_max=3)

DS_x_0([0.0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.5,2.0,2.5,3.0,3.5,4.0],"none",r_max=3) #[0.0,0.2,0.4,0.5,0.7,0.8,0.9,1.0,1.5,2.0,2.5,3.0,3.5,4.0],[0.0,0.2,0.5,1.1,1.25,1.3,1.4,1.5,1.6,2.0,2.5,4.0]
#D_nmax("bcc",f=0.7)
#sus(0.2, collect(0.05:0.01:3.0))

#DS_Bonca(collect(0.1:0.01:1.0),r_max=3) #only for square lattice

#wlim = 0.001
#xlim = 0.01
#DS_x_a(collect(0.01:0.01:1.0),wlim,[0.0,0.2,0.4,0.5],r_max=3)



#=
# check moments
c_iipDyn_mat_subst = load_ciipDyn_mat_subst(0.0)
moments = get_moments_Gklim(c_iipDyn_mat_subst, hte_lattice, ex)
c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat_subst)
Tχ = Polynomial([(-1)^n*sum(c_iipEqualTime_mat[i,1][n+1] for i in 1:hte_lattice.lattice.length) for n in 0:n_max])
#println(moments)
=#
#=
# u-Pade
x_vec_plot = collect(0.05:0.01:3.0)
f = 0.3
ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
u_vec = tanh.(f .* x_vec_plot)

plt = plot(plt, 1 ./ x_vec_plot, Tχ.(x_vec_plot))
plot!(plt, 1 ./ x_vec_plot, get_upade(Tχ,f=f,m=5,n=4).(u_vec))
#scatter!(plt,[0.5,1.0,1.5,2.0],get_upade(Tχ,f=f,m=5,n=4).(tanh.(f.*[0.5,1.0,1.5,2.0])))
display(plt)
#println(Tχ)
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