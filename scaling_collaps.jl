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
lattice_type = "triang"

# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = false
j3 = false
j4 = false

L = 12

hte_lattice = getLattice(L,lattice_type,j1,j2,j3,j4);

Γ,K,M = (0,0), (2*π/3,2*π/sqrt(3)), (0,2*π/sqrt(3))

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

###### dynamical Matsubara correlator (k-space)
k,k_label = K,"K"
a = 0.0 #0.06
sc = sqrt(1+a^2)
c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a)
c_kDyn_mat = get_c_k([k],c_iipDyn_mat_subst,hte_lattice)[1]
m_vec = get_moments_from_c_kDyn(c_kDyn_mat)

x_vec_bare = collect(0.0:0.025:1.2) # J1/T
x_vec = collect(0.0:0.2:4.0) #J1/T

x0_vec = 1 ./ [3.0,1.8,1.2,0.95,0.8,0.7,0.6,0.5,0.43,0.38]  # for these T/J1 the DSF will be computed

##### plot DSF and related quantities
if true
    w_vec = collect(0.0:0.02:3.7) # ω/J1

    ###### Pade for moments with x-series and u-series
    r_max = 3   # maximal order for moment
    r_min = 3
    intercept0 = true
    f = 0.55 # 0.55,0.6
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max+1,f)
    u_vec = tanh.(f .* x_vec)
    u0_vec = tanh.(f .* x0_vec)
    m0_vec = [Float64[] for _ in x0_vec]

    plt_m = Plots.plot([0],[0],xlims=(0,x_vec[end]),label="")
    Plots.plot!(plt_m,xlabel=L"x=J/T",ylabel=L"x \cdot m_{\mathbf{k},2r}(x) \, / \, m_{\mathbf{k},2r}(0)",legend=:topleft)
    Plots.plot!(plt_m,-x_vec,0*x_vec,color=:grey,label="x bare")
    Plots.plot!(plt_m,-x_vec,0*x_vec,color=:grey,linestyle=linestyle_vec[2],label="u Padé [7-r,6-r]")
    Plots.plot!(plt_m,-x_vec,0*x_vec,color=:grey,linestyle=linestyle_vec[3],label="u Padé [6-r,5-r]")
    Plots.annotate!(plt_m,3,2,Plots.text(L"\mathbf{k}="*string(k_label)*",  f="*string(f),7))
    
    for r in 0:r_max
        println("r=$r")
        poly_x = Polynomial([0,1],:x)
        xm_norm_r = coeffs(poly_x * (m_vec[1+r]/m_vec[1+r](0)))
        p_x = Polynomial(xm_norm_r)
        p_u = Polynomial(ufromx_mat[1:n_max+2-2*r,1:n_max+2-2*r]*xm_norm_r)
        
        Plots.plot!(plt_m,x_vec_bare.*sc,Polynomial(xm_norm_r).(x_vec_bare),color=color_vec[r+1],linewidth=0.4,label="r=$r",alpha=0.7) 

        ### x-Padé moments are not well behaved (grow large or negative)
        #Plots.plot!(plt_m,x_vec,get_pade(p_x,7-r,6-r).(x_vec),color=color_vec[r+1],linestyle=linestyle_vec[2],label="",alpha=0.7)
        #Plots.plot!(plt_m,x_vec,get_pade(p_x,6-r,5-r).(x_vec),color=color_vec[r+1],linestyle=linestyle_vec[3],label="",alpha=0.7)

        Plots.plot!(plt_m,x_vec.*sc,get_pade(p_u,7-r,6-r).(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[2],label="",alpha=0.7)
        Plots.plot!(plt_m,x_vec.*sc,get_pade(p_u,6-r,5-r).(u_vec),color=color_vec[r+1],linestyle=linestyle_vec[3],label="",alpha=0.7)

        ### extract moments from u-Padé [7-r,6-r] at x0_vec 
        for x0_pos in eachindex(x0_vec)
            x0 = x0_vec[x0_pos]
            u0 = u0_vec[x0_pos]
            push!(m0_vec[x0_pos],m_vec[1+r](0)/x0 * get_pade(p_u,7-r,6-r)(u0))
        end
    end

    ###### δ_r, JS and A for x ∈ x0_vec
    plt_δ = Plots.plot([0],[0],label="",xlabel=L"r",ylabel=L"\delta_{\mathbf{k},r}",legend=:bottomright)
    plt_JS = Plots.plot(xlims=(0,w_vec[end]),xlabel=L"\omega/J=w",ylabel=L"J \, S(\mathbf{k}="*k_label*L",\omega)",legendfontsize=5.0,legend=:topright)
    plt_JAo2π = Plots.plot(xlabel=L"\omega/J=w",ylabel=L"J \, A(\mathbf{k}="*k_label*L",\omega)",legend=:topleft)

    #plt_JSw0=Plots.plot([0.55,0.55],[0.0,0.09],xscale=:log10,xlims=(0.1,1.02/x0_vec[1]),label="roton-like energy [Zheng2006]",color=:grey,xlabel=L"T/J=1/x",ylabel=L"J \, S(\mathbf{k}="*k_label*L",\omega \rightarrow 0)",legend=:bottomright)

    for x0_pos in eachindex(x0_vec)
        x0 = x0_vec[x0_pos]

        ### plot Dyn-HTE
        δ_vec,r_vec = fromMomentsToδ(m0_vec[x0_pos])
        Plots.scatter!(plt_δ,r_vec,δ_vec,color=thermalCol13_vec[x0_pos],label="")
        
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

        δ_vec_ext = extrapolate_δvec(δ_vec,r_max,r_max,4000,true)
        Plots.plot!(plt_δ,r_max+1:6,δ_vec_ext[r_max+2:7],label="",color=thermalCol13_vec[x0_pos])

        JSw_vec = [JSwithTerminator(δ_vec,1.0*x0,w,extrap_params) for w in w_vec].*sc

        if k_label=="K"
            Plots.plot!(plt_JS,w_vec./sc,JSw_vec,color=thermalCol13_vec[x0_pos],label="")
        else
            Plots.plot!(plt_JS,w_vec./sc,JSw_vec,color=thermalCol13_vec[x0_pos],label="T/J="*string(round(1/x0,digits=2)))
        end

        Plots.plot!(plt_JAo2π,w_vec./sc, JSw_vec .* (1 .- exp.(-x0 .* w_vec)),color=thermalCol13_vec[x0_pos],label="")
        #Plots.scatter!(plt_JSw0,[1/x0],[JSw_vec[1]],color=thermalCol15_vec[x0_pos],label="")
    end

    ### plot results
    xPlots,yPlots=3,1
    plt_final = Plots.plot(plt_m, plt_δ, plt_JS, layout=(yPlots,xPlots), size=(0.7*aps_width*xPlots,0.7*aps_width*yPlots))
    display(plt_final)
    #Plots.savefig(plt_final,"CaseStudy/Triangular_Lattice/Triangular_DSF_k"*k_label*".png")
end

if false
    ### prepare insets
    if k_label=="K" 
        Plots.plot!(plt_JS,inset=Plots.bbox(0.62,0.03,0.36,0.49),subplot=2) 
    end
    if k_label=="M" 
        Plots.plot!(plt_JS,xlabel="",xformatter=_->"") 
        Plots.plot!(plt_JS,inset=Plots.bbox(0.2,0.65,0.3,0.3),subplot=2)
        Plots.plot!(plt_JS[2],xlims=(-1.3,1.3),ylims=(-1.2,1.2),aspect_ratio = :equal,xaxis=false,yaxis=false)
        a1 = [1/2, sqrt(3)/2]
        a2 = [1, 0]
        a3 = [-1/2, sqrt(3)/2]
        for x in -5:1:5, y in -5:1:5
            r = x*a1 .+ y*a2
            r1 = r .+ a1
            r2 = r .+ a2
            r3 = r .+ a3
            Plots.plot!(plt_JS[2],[r[1],r1[1]],[r[2],r1[2]],markers=:dot,color=:black,label="")
            Plots.plot!(plt_JS[2],[r[1],r2[1]],[r[2],r2[2]],markers=:dot,color=:black,label="")
            Plots.plot!(plt_JS[2],[r[1],r3[1]],[r[2],r3[2]],markers=:dot,color=:black,label="")
        end
    end

    ### scaling plot of DSF at k=K (as inset)
    if k_label=="K" && true
        w_max = 1.0
        w_vec = collect(0.0:0.02:w_max)
        α = 1.1

        Plots.annotate!(plt_JS[2],1,0.1,Plots.text(L"\alpha="*string(α)*"0(2)",7))
        Plots.plot!(plt_JS[2],xlims=(0,2.7),ylims=(0.08,0.165), guidefont=(6,), tickfont=(6,),xlabel=L"\omega/T\;\;(\omega \leq"*string(w_max)*L" \,J)",ylabel=L"J\, S(\mathbf{k}="*k_label*L",\omega) \cdot (T/J)^\alpha")

        for x0_pos in eachindex(x0_vec)
            x0 = x0_vec[x0_pos]
            if 1/x0 >= 0.3 && 1/x0 <= 1.0

                ### plot Dyn-HTE
                δ_vec,r_vec = fromMomentsToδ(m0_vec[x0_pos])
                Plots.scatter!(plt_δ,r_vec,δ_vec,color=thermalCol13_vec[x0_pos],label="x=$x0")
                
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

                JSw_vec = [JSwithTerminator(δ_vec,1.0*x0,w,extrap_params) for w in w_vec].*sc

                Plots.plot!(plt_JS[2],w_vec .* (x0/sc), JSw_vec .* (1/x0)^α ,color=thermalCol13_vec[x0_pos],label="")
            end
        end
        display(plt_JS)
    end


    ### final plot for paper
    if k_label=="K" plt_JS_K = deepcopy(plt_JS) end
    if k_label=="M" plt_JS_M = deepcopy(plt_JS) end


    ###### run the above for both K and M and then put together
    xPlots,yPlots=1,2
    plt_final = Plots.plot(plt_JS_M,plt_JS_K, layout=(yPlots,xPlots), size=(aps_width*xPlots,(0.45)*aps_width*yPlots),dpi=600)
    display(plt_final)
    #Plots.savefig(plt_final,"CaseStudy/Triangular_Lattice/Triangular_DSF.png")
end

if k_label=="K" && true
    w_max = 1.0
    w_vec = collect(0.0:0.02:w_max)
    α1,α2,α3 = [1.0,1.1,1.2]

    plt_JS_scaled1 = Plots.plot(ylims=(0.08,0.27),legend=:topright)
    plt_JS_scaled2 = Plots.plot(ylims=(0.08,0.27),legend=:topright)
    plt_JS_scaled3 = Plots.plot(ylims=(0.08,0.27),legend=:topright)
    Plots.plot!(plt_JS_scaled1,xlabel=L"\omega/T\;\;(\omega \leq"*string(w_max)*L" \,J)",ylabel=L"J\, S(\mathbf{k}="*k_label*L",\omega) \times (T/J)^\alpha  \;\;\; \alpha="*string(α1))
    Plots.plot!(plt_JS_scaled2,xlabel=L"\omega/T\;\;(\omega \leq"*string(w_max)*L" \,J)",ylabel=L"J\, S(\mathbf{k}="*k_label*L",\omega) \times (T/J)^\alpha  \;\;\; \alpha="*string(α2))
    Plots.plot!(plt_JS_scaled3,xlabel=L"\omega/T\;\;(\omega \leq"*string(w_max)*L" \,J)",ylabel=L"J\, S(\mathbf{k}="*k_label*L",\omega) \times (T/J)^\alpha  \;\;\; \alpha="*string(α3))


    for x0_pos in eachindex(x0_vec)
        x0 = x0_vec[x0_pos]
        if 1/x0 >= 0.3 && 1/x0 <= 1.0

            ### plot Dyn-HTE
            δ_vec,r_vec = fromMomentsToδ(m0_vec[x0_pos])
            Plots.scatter!(plt_δ,r_vec,δ_vec,color=thermalCol13_vec[x0_pos],label="x=$x0")
            
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

            JSw_vec = [JSwithTerminator(δ_vec,1.0*x0,w,extrap_params) for w in w_vec].*sc

            Plots.plot!(plt_JS_scaled1,w_vec .* (x0/sc), JSw_vec .* (1/x0)^α1 ,color=thermalCol13_vec[x0_pos],label="T/J="*string(round(1/x0,digits=2)))
            Plots.plot!(plt_JS_scaled2,w_vec .* (x0/sc), JSw_vec .* (1/x0)^α2 ,color=thermalCol13_vec[x0_pos],label="T/J="*string(round(1/x0,digits=2)))
            Plots.plot!(plt_JS_scaled3,w_vec .* (x0/sc), JSw_vec .* (1/x0)^α3 ,color=thermalCol13_vec[x0_pos],label="T/J="*string(round(1/x0,digits=2)))

        end
    end

    xPlots,yPlots=3,1
    plt_final = Plots.plot(plt_JS_scaled1,plt_JS_scaled2,plt_JS_scaled3, layout=(yPlots,xPlots), size=(aps_width*xPlots,0.5*aps_width*yPlots))
    display(plt_final)
    #Plots.savefig(plt_final,"CaseStudy/Triangular_Lattice/Triangular_JS_k"*k_label*"_scaling.png")
end