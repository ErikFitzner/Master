using JLD2, DelimitedFiles, SimpleWeightedGraphs, Plots, Symbolics
using Profile, ProfileView
using FlameGraphs, FileIO

include("plotConventions.jl")
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 

### load graph evaluations
spin_length = 1/2
n_max = 10

### prepare lattice
lattice_type = "Shastry-Sutherland"

# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = true
j3 = false
j4 = false

L = 10

hte_lattice = getLattice(L,lattice_type,j1,j2,j3,j4);

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
    fileName_c = "CaseStudy/$(lattice_type)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_J1_$(1*j1)_J2_$(1*j2)_J3_$(1*j3)_J4_$(1*j4).jld2"
    if isfile(fileName_c)
        println("loading "*fileName_c)
        c_iipDyn_mat = load_object(fileName_c)
    else
        hte_graphs, C_Dict_vec = getGraphsG(spin_length,n_max)
        #Profile.clear()
        #=@profile=# c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs,C_Dict_vec);
        #ProfileView.view()
        save_object(fileName_c,c_iipDyn_mat)
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


###### load an existing flamegraph
#=
data, lidict = load("C:/Users/User/Downloads/Profile_cii_2.jlprof")   # tuple zurück
g = flamegraph(data; lidict)                   # Flame-Graph bauen
ProfileView.view(g)                            # interaktiv inspizieren
=#
