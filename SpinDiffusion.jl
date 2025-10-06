using JLD2, DelimitedFiles, SimpleWeightedGraphs, Plots, Symbolics

include("plotConventions.jl")
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 

@variables x1 x2 Δ

### load graph evaluations
spin_length = 1/2
bc1 = spin_length*(spin_length+1)/3
n_max = 6

### prepare lattice
lattice_type = "simple_cubic" #"Shastry-Sutherland"
d = 3
ex = zeros(d)  # creates a d-dimensional zero vector
ex[1] = 1      # set first component to 1

# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = true
j3 = false
j4 = false

L = 6

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

function Tchi(c_iipDyn_mat)
    c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat)
    result_suscept = sum([get_TGiip_Matsubara_xpoly(c_iipEqualTime_mat, j, 1, 0) for j in 1:hte_lattice.lattice.length])
    return substitute(result_suscept, Dict(x2 => 0))
end

function compute_Gklim(c_iipDyn_mat, hte_lattice, ex)
    lattice = hte_lattice.lattice
    center_sites = hte_lattice.basis_positions
    z = 0
    for b in 1:length(center_sites)
        for i in 1:length(lattice)
            z += (dot(ex, getSitePosition(lattice,i) .- getSitePosition(lattice,center_sites[b])))^2 * get_TGiip_Matsubara_xpoly(c_iipDyn_mat, i, b, 1)
        end
    end
    return (-1//2) * z / length(center_sites)
end

function extract_md(expr)
    expr_expanded = SymbolicUtils.expand(expr)
    terms = Symbolics.istree(expr_expanded) && Symbolics.operation(expr_expanded) == (+) ? Symbolics.arguments(expr_expanded) : (expr_expanded,)
    md_dict = Dict{Int, Any}()
    for term in terms
        deg_Δ = Symbolics.degree(term, Δ)
        if iseven(deg_Δ)
            r = div(deg_Δ, 2)
            coeff_val = simplify(-(-1)^r * term / (x1^(2r) * Δ^(2r)))
            md_dict[r] = get(md_dict, r, 0) + coeff_val
        end
    end
    return md_dict
end

function Diffusion_Gauß(Tχ, Gklim_dict, x)
    m2 = Symbolics.build_function(Gklim_dict[1], x1; expression=Val{false})(x)
    m4 = Symbolics.build_function(Gklim_dict[2], x1; expression=Val{false})(x)
    Tχ = Symbolics.build_function(Tχ, x1; expression=Val{false})(x)
    result = sqrt(pi)/(Tχ*sqrt(2))*sqrt(m2^3/m4)

    return result
end

function calc_diff(c_iipDyn_mat, hte_lattice, ex, x_vec)
    Gklim = compute_Gklim(c_iipDyn_mat, hte_lattice, ex)
    Gklim_subst = Symbolics.value(Symbolics.expand(substitute(Gklim, Dict(x2 => 0))))
    Gklim_dict = extract_md(Gklim_subst)
    Tχ = Tchi(c_iipDyn_mat)
    return [Diffusion_Gauß(Tχ,Gklim_dict,x) for x in x_vec]
end

# plot
if true
    x_vec = collect(0.0:0.01:1.0)
    y_vec = calc_diff(c_iipDyn_mat, hte_lattice, ex, x_vec)
    plt = Plots.plot([0],[0],label="",xlabel=L"J/T",ylabel=L"D_S/J",title="Spin Diffusion, $(lattice_type)",size=(1.5*aps_width,aps_width), legend=:topleft)
    Plots.plot!(plt, x_vec, y_vec, label="Dyn-HTE")
    
    #### hypercubic lattice
    #Plots.scatter!(plt, [0], [sqrt(pi*bc1/(4*d-2-1/(4*bc1)))], label="[Kopietz1993]", markersize=6)
    
    display(plt)
    #savefig(plt,"Images/SpinDiffusion_$(lattice_type).png")
end
