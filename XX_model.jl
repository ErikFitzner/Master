using JLD2, DelimitedFiles, SimpleWeightedGraphs, Plots, Symbolics

#include("plotConventions.jl")
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 

### load graph evaluations
spin_length = 1/2
n_max = 4

### prepare lattice
lattice_type = "square"

# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = true
j3 = false
j4 = false

L = 4

hte_lattice = getLattice(L,lattice_type,j1,j2,j3,j4);
#println(ne(hte_lattice.graph))
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
if true
    fileName_c = "CaseStudy/$(lattice_type)_XX_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_J1_$(1*j1)_J2_$(1*j2)_J3_$(1*j3)_J4_$(1*j4).jld2"
    if isfile(fileName_c)
        println("loading "*fileName_c)
        c_iipDyn_mat = load_object(fileName_c)
    else
        hte_graphs, C_Dict_vec = getGraphsG_XX(spin_length,n_max)
        #Profile.clear()
        #=@profile=# c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs,C_Dict_vec, verbose=true);
        #ProfileView.view()
        save_object(fileName_c,c_iipDyn_mat)
    end
end

if true
    a_vec = [0.0] #[0.0,0.08,0.2,0.35,0.4,0.45,0.5,0.7,1.0] # [0.0,0.2,0.4,0.5,0.7,0.8,0.9,1.0,1.5,2.0,2.5,3.0,3.5,4.0] # [0.0,0.047,0.08,0.12,0.25,0.5,1.0,2.0]
    b = 0.0
    c = 0.0

    for a in a_vec
        fileName = "CaseStudy/$(lattice_type)_XX_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_a_$(a)_b_$(b)_c_$(c).jld2"
        println("substituting c_iipDyn_mat with a=$(a), b=$(b), c=$(c)")
        c_iipDyn_mat_subst = get_c_iipDyn_mat_subst(c_iipDyn_mat,hte_lattice,a,b,c);
        save_object(fileName,c_iipDyn_mat_subst)
        if false
            fileName2 = "CaseStudy/$(lattice_type)_XX_collapse_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_a_$(a)_b_$(b)_c_$(c).jld2"
            println("Collapsing ladder to the chain")
            c_iipDyn_mat_subst_collaps = Array{Matrix{Float64}}(undef, 2*L+1, 1)
            for i in 1:2*L+1
                j = 2*i-1
                c_iipDyn_mat_subst_collaps[i,1] = 1/4 * (c_iipDyn_mat_subst[j,1]+c_iipDyn_mat_subst[j,2]+c_iipDyn_mat_subst[j+1,1]+c_iipDyn_mat_subst[j+1,2])
            end
            save_object(fileName,c_iipDyn_mat_subst_collaps)
        end
    end
end

###### check results
if false
    # Matsubara
    #result = get_TGiip_Matsubara_xpoly(c_iipDyn_mat,36,1,0)
    #println("result = $(result)")

    # Uniform susceptebility
    @variables x1 x2
    c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat)
    result_suscept = sum([get_TGiip_Matsubara_xpoly(c_iipEqualTime_mat, j, 1, 0) for j in 1:hte_lattice.lattice.length])
    println("uniform susceptebility = $(result_suscept)")
    expr_sub = substitute(result_suscept, Dict(x2 => 0))
    println("uniform susceptibility (x2=0) = $(expr_sub)")
end
