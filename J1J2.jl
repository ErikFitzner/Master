using JLD2, DelimitedFiles, SimpleWeightedGraphs, Plots, Symbolics
#using Profile, ProfileView
include("plotConventions.jl")
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 

# TODO add an if statement to distinguish case using one bond type from case with different types --> unique_gG_vec vs GraphG
#@variables x1 x2
### load graph evaluations
spin_length = 1/2
n_max = 9

### prepare lattice
lattice_type = "Shastry-Sutherland"
#Γ,K,M = (0,0), (2*π/3,2*π/sqrt(3)), (0,2*π/sqrt(3))

# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = true
j3 = false
j4 = false

L = 9

hte_lattice = getLattice(L,lattice_type,j1,j2,j3,j4);
#println(hte_lattice.basis_positions)

### plot lattice
if false
    weights = hte_lattice.graph.weights
    # Assign a color for each bond weight: 1=blue, 2=red, 3=green, 4=orange, else=gray
    color_map = Dict(1.0 => :blue, 2.0 => :red, 3.0 => :green, 4.0 => :orange)
    edge_colors = [get(color_map, w, :gray) for w in weights]
    display(graphplot(hte_lattice.graph, names=1:nv(hte_lattice.graph), markersize=0.2, fontsize=8, nodeshape=:rect, curves=false, edgecolor=edge_colors))
end

### compute all correlations in the lattice (or load them)
if false
    fileName_c = "CaseStudy/$(lattice_type)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_J1_$(1*j1)_J2_$(1*j2)_J3_$(1*j3)_J4_$(1*j4).jld2"
    if isfile(fileName_c)
        println("loading "*fileName_c)
        c_iipDyn_mat = load_object(fileName_c)
    else
        hte_graphs, C_Dict_vec = getGraphsG(spin_length,n_max)
        c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs,C_Dict_vec);
        save_object(fileName_c,c_iipDyn_mat)
    end
end

#test = load_object("CaseStudy/Triangular_Lattice/Triangular_Lattice_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(12) * "_L" * string(12) * ".jld2")
#println("test = $(test[1,1])")

#println(size(c_iipDyn_mat))

#result = get_TGiip_Matsubara_xpoly(c_iipDyn_mat,36,1,0)
#println("result = $(result)")

#c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat)
#result_suscept = sum([get_TGiip_Matsubara_xpoly(c_iipEqualTime_mat, j, 1, 0) for j in 1:hte_lattice.lattice.length])
#println("uniform susceptebility = $(result_suscept)")
#expr_sub = substitute(result_suscept, Dict(x2 => 0))
#println("uniform susceptibility (x2=0) = $(expr_sub)")

if false
a = 0.5 #J/J_D
f_num = subvalue(result_suscept,1/a)
xs1 = range(0.1, 2, length=200)
xs2 = range(0.1, 10, length=200)
ys = [f_num(x) for x in xs1]
y_pade = robustpade(f_num,6,6).(xs2)

Plots.plot(1 ./ xs1, ys, xlabel=L"T/J", ylabel="χ", title=L"J/J_D="*string(a), label="x-Series")
#Plots.plot!(1 ./ xs2, y_pade, color="orange", linestyle=:dash, alpha=0.7, label="Pade")
display(current())
end

#########################################################################################
###### Dynamic structure factor (DSF) ###################################################
#########################################################################################
if false
    #k_pihalf = []
    #a_vec = 1 ./ [0.5,0.8,1.0]
    #f_vec = [0.7,0.5,0.5]
    #substitute specific values for x2,x3,x4 to reduce to the previous form of c_iipDyn_mat
    #for i in 1:length(a_vec)
    a = 0.0 #1/0.5 # J_2/J_1
    #a = a_vec[i]
    b = 0.0
    c = 0.0

    fileName_c = "CaseStudy/$(lattice_type)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_a_$(a)_b_$(b)_c_$(c).jld2"
    if isfile(fileName_c)
        println("loading "*fileName_c)
        c_iipDyn_mat_subst = load_object(fileName_c)
    else
        println("substituting c_iipDyn_mat with a=$(a), b=$(b), c=$(c)")
        c_iipDyn_mat_subst = get_c_iipDyn_mat_subst(c_iipDyn_mat,hte_lattice,a,b,c);
        save_object(fileName_c,c_iipDyn_mat_subst)
    end

    w_vec = collect(0.0:0.025:3.5)
    r_max = 3                
    f = 0.7
    #f = f_vec[i]
    ufromx_mat = get_LinearTrafoToCoeffs_u(n_max+1,f)
    poly_x = Polynomial([0,1],:x)

    x0 = 2.0
    u0 = tanh.(f .* x0)

    ### define and generate k-path 
    path = [(pi/2,pi/2),(pi,0),(pi,pi),(pi/2,pi/2),(0.001,0.001),(pi,0)]
    pathticks = ["(π/2,π/2)","(π,0)","(π,π)","(π/2,π/2)","(0,0)","(π,0)"]

    Nk = 75  #75
    k_vec,kticks_positioins = create_brillouin_zone_path(path, Nk)
    JSkw_mat = zeros(Nk+1,length(w_vec))

    ### fill JSkw_mat
    for k_pos in eachindex(k_vec)
        k = k_vec[k_pos]

        c_kDyn = get_c_k(k,c_iipDyn_mat_subst,hte_lattice)
        m_vec = get_moments_from_c_kDyn(c_kDyn)
        m0 = Float64[]

        for r in 0:r_max
            xm_norm_r = coeffs(poly_x * (m_vec[1+r]/m_vec[1+r](0)))
            p_u = Polynomial(ufromx_mat[1:n_max+2-2*r,1:n_max+2-2*r]*xm_norm_r)
            
            push!(m0,m_vec[1+r](0)/x0 * get_pade(p_u,7-r,6-r)(u0))
        end

        δ_vec,r_vec = fromMomentsToδ(m0)
        δ_vec_ext = extrapolate_δvec(δ_vec,length(δ_vec)-1,length(δ_vec)-1,2000,true)
        JSkw_mat[k_pos,:] = [JS(δ_vec_ext,1.0*x0,w,0.02) for w in w_vec]
    end
    #push!(k_pihalf,JSkw_mat[1,:])
    #end
end

if false
    a_vec = 1 ./ [0.5,0.8,1.0]
    w_vec = collect(0.0:0.025:3.5)
    plt = Plots.plot(xlabel=L"\omega/J_1",ylabel=L"J_1S(k,\omega)",legend=:topleft)
    for i in 1:length(k_pihalf)
    Plots.plot!(plt,w_vec,k_pihalf[i],label=L"J_2/J_1="*string(a_vec[i])*", "*L"k=(\pi/2,\pi/2)")
    end
    display(plt)
    #savefig(plt,"Images/Shastry-Sutherland_(pihalf,pihalf).png")
end

### plot JS(k,ω)
if true
    using CairoMakie

    fig = Figure(fontsize=8,size=(aps_width,0.6*aps_width));
    ax=Axis(fig[1,1],xlabel=L"\mathbf{k}",ylabel=L"\omega/J_1=w",xlabelsize=8,ylabelsize=8);
    hm=CairoMakie.heatmap!(ax,collect(0:Nk)/(Nk),w_vec, JSkw_mat,colormap=:viridis,colorrange=(0.001,0.4),highclip=:white);
    ax.xticks = ((kticks_positioins .- 1)/(Nk),pathticks)
    CairoMakie.Colorbar(fig[:, end+1], hm,size=8, label = L"J_1 S(\mathbf{k},\omega)")
    CairoMakie.text!(ax,"x=J1/T="*string(x0),position=[(0.05,0.8)],color=:white)
    CairoMakie.text!(ax,"J2/J1="*string(a),position=[(0.05,0.5)],color=:white)
    CairoMakie.text!(ax,"f=$f",position=[(0.05,0.2)],color=:white)

    resize_to_layout!(fig);
    display(fig)

    #save("Images/Shastry_Sutherland_JSkw_a_$(a)_2.png",fig; px_per_unit=6.0)
end