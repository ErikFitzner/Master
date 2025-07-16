using JLD2, DelimitedFiles, SimpleWeightedGraphs, Plots, Symbolics
include("plotConventions.jl")
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 

# TODO add an if statement to distinguish case using one bond type from case with different types --> unique_gG_vec vs GraphG 
# TODO get_c_iipDyn_mat mapreduce wieder einbauen?

### load graph evaluations
spin_length = 1/2
n_max = 4

### prepare lattice
lattice_type = "chain"
#Γ,K,M = (0,0), (2*π/3,2*π/sqrt(3)), (0,2*π/sqrt(3))

# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = true
j3 = false
j4 = false

L = 4

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
fileName_c = "CaseStudy/$(lattice_type)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_J1_$(1*j1)_J2_$(1*j2)_J3_$(1*j3)_J4_$(1*j4).jld2"
if isfile(fileName_c)
    println("loading "*fileName_c)
    c_iipDyn_mat = load_object(fileName_c)
else
    hte_graphs, C_Dict_vec = getGraphsG(spin_length,n_max)
    c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs,C_Dict_vec);
    #c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice.graph, hte_lattice.basis_positions, hte_graphs, C_Dict_vec);
    save_object(fileName_c,c_iipDyn_mat)
end

#test = load_object("CaseStudy/Triangular_Lattice/Triangular_Lattice_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(12) * "_L" * string(12) * ".jld2")
#println("test = $(test[1,1])")

#result = get_TGiip_Matsubara_xpoly(c_iipDyn_mat,1,1,0)
#println("result = $(result)")

#c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat)
#result_suscept = sum([get_TGiip_Matsubara_xpoly(c_iipEqualTime_mat, j, 1, 0) for j in 1:length(c_iipDyn_mat)])
#println("uniform susceptebility = $(result_suscept)")

# Substitute specific values
if false
    m = 1
    a = 0.6
    @variables x1 x2 x3 x4 Δ
    vars_in_result = Symbolics.get_variables(result)
    subs_dict = Dict(
        #x1 => 4,
        x2 => a*x1,
        #x3 => x1,
        #x4 => 0,
        #Δ => 1/(2*pi*m)
    )
    filtered_subs = Dict(k => v for (k,v) in subs_dict if any(isequal(k), vars_in_result))
    subs_expr = substitute(result, filtered_subs)
    # Simplify or evaluate numerically
    subs_expr_val = Symbolics.value(subs_expr)
    #println("specific values: $(subs_expr_val)")

    # Check if x1 is in the expression
    if occursin("x1", string(subs_expr_val))
        # Build numerical function f_num(x1)
        f_num = Symbolics.build_function(subs_expr_val, x1; expression=Val{false})
        xs1 = range(0, 2.5, length=250)
        xs2 = range(0, 4, length=400)
        ys = [f_num(x) for x in xs1]
        y_pade = robustpade(f_num,6,6).(xs2)
        y_pade2 = robustpade(f_num,5,5).(xs2)
        y_exact = [exact11(x, a) for x in xs2]

        # Prepare u-Padé
        g = 0.35
        u_vec = tanh.(g .* xs2)
        # Define f_num_u(u) = f_num(atanh(u)/g)
        f_num_u(u) = f_num(atanh(u)/g)
        y_pade_u = robustpade(f_num_u,6,6).(u_vec)
        y_pade_u2 = robustpade(f_num_u,5,5).(u_vec)
        
        # Plot
        plot(xs1, ys, xlabel="x1", ylabel="TG_11(m=0)", title="x2=$(a)*x1")
        plot!(xs2, y_pade, color="red", alpha=0.7, label="x-Padé[6,6]", linestyle=:dash)
        #plot!(xs2, y_pade2, color="red", alpha=0.7, label="x-Padé[5,5]", linestyle=:dashdot)
        plot!(xs2, y_exact, color="green", alpha=0.7, label="exact")
        
        plot!(xs2, y_pade_u, color="orange",linestyle=:dash,alpha=0.7,label="u-Padé[6,6] (f=$g)")
        plot!(xs2, y_pade_u2, color="orange",linestyle=:dashdot,alpha=0.7,label="u-Padé[5,5] (f=$g)")
        
        display(current())
        savefig("Images/G_11_J2=$(a)J1_u.png")
    else
        println("subs_expr_val does not depend on x1, skipping plot.")
    end
end

#########################################################################################
###### Dynamic structure factor (DSF): k-path through BZ ############################
#########################################################################################

###SPIN STRUCTURE FACTOR HEATMAPS
using CairoMakie
x = 2.0                #define temperature (x=J/T)
k_step_size = 1/41     #define k step size (in 1/π)
w_step_size = 0.025    #define ω step size (in 1/J)
#define k and ω vectors 
k_vec = vcat(vcat((0.0001,0.0),[(k*pi,0.0) for k in 0:k_step_size:2][2:end-1] ),(1.999*pi,0.0))#[(k,0.0) for k in 0.01:0.0039*2.4:(2*π-0.01)]
w_vec = collect(-3:w_step_size:3)

#calculate the spin structure factor for the given k and ω 
JSkw_mat = get_JSkw_mat("u_pade",x,k_vec,w_vec,c_iipDyn_mat,hte_lattice,r_min=3,r_max=3,r_ext=1000,f=0.48)

#plot the result
fig = Figure(size=(400,400),fontsize=20)
ax=Axis(fig[1,1],limits=(0,2,-3,3),xlabel=L"k/ \pi",ylabel=L"\omega/J=w",title="x="*string(x),titlesize=20,xlabelsize=20,ylabelsize=20)
hm=CairoMakie.heatmap!(ax,[k[1]/π for k in k_vec],w_vec,JSkw_mat,colormap=:viridis,colorrange=(0.0,0.45),highclip=:white)

#SHOW PLOT
##############
display(fig)