using JLD2, DelimitedFiles, SimpleWeightedGraphs, Plots, Symbolics, Richardson

include("plotConventions.jl")
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 

@variables x1 x2 Δ

### load graph evaluations
spin_length = 1/2
bc1 = spin_length*(spin_length+1)/3
n_max = 12

### prepare lattice
lattice_type = "chain" #"Shastry-Sutherland"

# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = true
j3 = false
j4 = false

L = 12

hte_lattice = getLattice(L,lattice_type,j1,j2,j3,j4);

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

function get_QFI(moments)
    if length(moments) < 7
        moments = vcat(moments, zeros(Int, 7 - length(moments)))
    end
    x_poly = Polynomial([0, 1]) 
    f_Q_poly = 4 -
    x_poly^2/3  * moments[2] +
    x_poly^4/60 * moments[3] -
    x_poly^6/1512 * moments[4] +
    x_poly^8/43200 * moments[5] -
    x_poly^10/1330560 * moments[6] +
    691 * x_poly^12/29719872000 * moments[7]
    return pi * f_Q_poly
end

#calculate the moments
a = 0.0
c_iipDyn_mat_subst = load_ciipDyn_mat_subst(a)
k = (pi/2,0) # (4*pi/3,0)
c_kDyn_mat = get_c_k([k],c_iipDyn_mat_subst,hte_lattice)[1]
m_vec = get_moments_from_c_kDyn(c_kDyn_mat)
f_Q_poly = get_QFI(m_vec)

x_vec = collect(0.1:0.01:7)
f = 0.23
ufromx_mat = get_LinearTrafoToCoeffs_u(n_max,f)
u_vec = tanh.(f .* x_vec)
coeffs_x = coeffs(f_Q_poly)
p_u = Polynomial(ufromx_mat*coeffs_x)

plt = plot(xlabel="J/T", ylabel=L"I_S[k,T]", title="$(lattice_type), a=$(a)", legend=:bottomright)
Plots.plot!(plt,x_vec,robustpade(p_u,5,6).(u_vec),color="blue",linestyle=:dashdot,alpha=0.7,label="u-Padé[5,6] (f=$f)")
Plots.plot!(plt,x_vec,robustpade(p_u,4,5).(u_vec),color="red",linestyle=:dash,alpha=0.7,label="u-Padé[4,5] (f=$f)")
#display(plt)
#savefig(plt,"Images/QFI/$(lattice_type)_a_$(a).png")

if true
    inv_x_vec = 1.0 ./ x_vec
    
    plt2 = plot(xlabel="T/J", ylabel=L"f_Q[k,T]", xlabelfontsize = 14, xtickfontsize = 12, ylabelfontsize = 14, ytickfontsize = 12, legendfontsize = 12)
    
    #Plots.plot!(plt2, inv_x_vec, f_Q_poly.(x_vec),color="gray",label="x-Series")
    Plots.plot!(plt2,inv_x_vec,robustpade(p_u,5,6).(u_vec),color="blue",linestyle=:dash,alpha=0.7,label="u-Padé[5,6] (f=$f)")
    #Plots.plot!(plt2,inv_x_vec,robustpade(p_u,4,5).(u_vec),color="red",linestyle=:dashdot,alpha=0.7,label="u-Padé[4,5] (f=$f)")
    display(plt2)
end