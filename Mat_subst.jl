using JLD2, DelimitedFiles, SimpleWeightedGraphs, Plots, Symbolics
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 

### load graph evaluations
spin_length = 1/2
n_max = 12

### prepare lattice
lattice_type = "Shastry-Sutherland"

# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = true
j3 = false
j4 = false

L = 12

hte_lattice = getLattice(L,lattice_type,j1,j2,j3,j4);

a_vec = [0.0,0.2,0.5,1.1,1.25,1.3,1.4,1.5,1.6,2.0,2.5,4.0]
b = 0.0
c = 0.0

fileName = "CaseStudy/$(lattice_type)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_J1_$(1*j1)_J2_$(1*j2)_J3_$(1*j3)_J4_$(1*j4).jld2"
if isfile(fileName)
    println("loading "*fileName)
    c_iipDyn_mat = load_object(fileName)
end

for a in a_vec
    fileName_c = "CaseStudy/$(lattice_type)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_a_$(a)_b_$(b)_c_$(c).jld2"
    println("substituting c_iipDyn_mat with a=$(a), b=$(b), c=$(c)")
    c_iipDyn_mat_subst = get_c_iipDyn_mat_subst(c_iipDyn_mat,hte_lattice,a,b,c);
    save_object(fileName_c,c_iipDyn_mat_subst)
end