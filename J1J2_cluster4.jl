using JLD2, DelimitedFiles, SimpleWeightedGraphs, Plots, Symbolics
#using Profile, ProfileView
include("plotConventions.jl")
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 

# TODO add an if statement to distinguish case using one bond type from case with different types --> unique_gG_vec vs GraphG

### load graph evaluations
spin_length = 1/2
n_max = 3

### prepare lattice
lattice_type = "cluster4"

# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = true
j3 = false
j4 = false

L = 4

hte_lattice = getLattice(L,lattice_type,j1,j2,j3,j4);

### plot lattice
if true
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
#result_suscept = sum([get_TGiip_Matsubara_xpoly(c_iipEqualTime_mat, j, 1, 0) for j in 1:length(c_iipDyn_mat[:,1])]) # for b in 1:length(c_iipDyn_mat[1,:])])
#println("uniform susceptebility = $(result_suscept)")

# Substitute specific values
if false
    a = 0.5
    f_num = subvalue(result,a)
    xs1 = range(0, 4, length=400)
    xs2 = range(0, 4, length=400)
    ys = [f_num(x) for x in xs1]
    y_pade = robustpade(f_num,6,6).(xs2)
    y_pade2 = robustpade(f_num,5,5).(xs2)
    #y_exact = [exact12(x, a) for x in xs2]

    # Prepare u-Padé
    g = 0.35
    u_vec = tanh.(g .* xs2)
    # Define f_num_u(u) = f_num(atanh(u)/g)
    f_num_u(u) = f_num(atanh(u)/g)
    y_pade_u = robustpade(f_num_u,6,6).(u_vec)
    y_pade_u2 = robustpade(f_num_u,5,5).(u_vec)
        
    # Plot
    Plots.plot(xs1, ys, xlabel="x1", ylabel=L"TG_{12}(m=0)", title="x2=$(a)*x1", label="x-Series")
    Plots.plot!(xs2, y_pade, color="red", alpha=0.7, label="x-Padé[6,6]", linestyle=:dash)
    Plots.plot!(xs2, y_pade2, color="red", alpha=0.7, label="x-Padé[5,5]", linestyle=:dashdot)
    #Plots.plot!(xs2, y_exact, color="green", alpha=0.7, label="exact")
        
    Plots.plot!(xs2, y_pade_u, color="orange",linestyle=:dash,alpha=0.7,label="u-Padé[6,6] (f=$g)")
    Plots.plot!(xs2, y_pade_u2, color="orange",linestyle=:dashdot,alpha=0.7,label="u-Padé[5,5] (f=$g)")
        
    display(current())
    #savefig("Images/G_12_J2=$(a)J1_u.png")
end