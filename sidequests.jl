using Plots, Symbolics
include("plotConventions.jl")
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 

###### Plot the movement of the DSF Maxima
if false
    # Position of the maxima
    α_list=[0.0,0.2,0.45,0.7,1.0]
    pos_dyn=[[(1,0.13)],[(1,0.13)],[(0.5,1),(1.5,1)],[(0.65,0.73),(1.35,0.73)],[(0.65,0.67),(1.35,0.67)]]
    pos_vmc=[[(1,0)],[(1,0)],[(0.57,1),(1.43,1)],[(0.58,0.67),(1.42,0.67)],[(0.5,0.42),(1.5,0.42)]]
        
    # Assign a color per α
    colors = palette(:viridis, length(α_list))

    # Make the base plot (adjust to your heatmap or other data plot)
    plt = plot(title="DSF Maxima Positions", xlabel="k/π", ylabel="ω/J₁", legend=true, xlim=(0,2), ylim=(0,3.5))

    # Plot maxima as scatter points
    for (i, pos) in enumerate(pos_dyn)
        x_vals = [p[1] for p in pos]
        y_vals = [p[2] for p in pos]
        scatter!(plt, x_vals, y_vals, color=colors[i], label="α = $(α_list[i])", markersize=6)
    end
    for (i, pos) in enumerate(pos_vmc)
        x_vals = [p[1] for p in pos]
        y_vals = [p[2] for p in pos]
        scatter!(plt, x_vals, y_vals, color=colors[i], marker=:cross, markersize=6, label=false)
    end
    display(plt)
    #savefig("Images/DSF_pos_max.png")
end

###### various checks for c_iipDyn_mat
if true
    spin_length = 1/2
    n_max = 12

    lattice_type1 = "aniso_square" #"Shastry-Sutherland"
    lattice_type2 = "aniso_square"

    j1 = true
    j2 = true
    j3 = false
    j4 = false

    L = 12

    hte_lattice1 = getLattice(L,lattice_type1,j1,j2,j3,j4);
    hte_lattice2 = getLattice(L,lattice_type2,j1,j2,j3,j4);

    #println(length(hte_lattice1.lattice.sitePositions))
    #println(hte_lattice2.lattice.sitePositions[hte_lattice2.basis_positions[4]])
    #println(hte_lattice2.lattice.sitePositions[31])
    #println(hte_lattice2.basis_positions)
    
    #fileName1 = "CaseStudy/$(lattice_type1)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_J1_$(1*j1)_J2_$(1*j2)_J3_$(1*j3)_J4_$(1*j4).jld2"
    fileName1 = "CaseStudy/$(lattice_type1)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_a_$(0.0)_b_$(0.0)_c_$(0.0).jld2"
    #fileName2 = "CaseStudy/$(lattice_type2)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_J1_$(1*j1)_J2_$(1*j2)_J3_$(1*j3)_J4_$(1*j4).jld2"
    fileName2 = "CaseStudy/$(lattice_type2)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_a_$(1.0)_b_$(0.0)_c_$(0.0).jld2"
    
    c_iipDyn_mat1 = load_object(fileName1)
    c_iipDyn_mat2 = load_object(fileName2)

    #=
    for i in eachindex(c_iipDyn_mat1)
        A = c_iipDyn_mat1[i]
        B = c_iipDyn_mat2[i]

        if !isapprox(A, B; atol=1e-6) #!isequal(A, B)
            println("Mismatch at index $i")
            println("  A = ", A)
            println("  B = ", B)
        end
    end
    =#

    #println(isequal(c_iipDyn_mat1,c_iipDyn_mat2))
    #println(isapprox(c_iipDyn_mat1, c_iipDyn_mat2; atol=1e-6))
    
    #println(size(c_iipDyn_mat1))
    #println(size(c_iipDyn_mat2))

    #println(c_iipDyn_mat1[12,1])
    #println(c_iipDyn_mat2[22,1])

    #println(get_c_iipDyn_mat_subst_old(c_iipDyn_mat1,hte_lattice1,0.3,0.0,0.0)[41])
    #println(get_c_iipDyn_mat_subst(c_iipDyn_mat2,hte_lattice2,0.0,0.0,0.0)[41])

    #=
    path = [(pi/2,pi/2),(pi,0),(pi,pi),(pi/2,pi/2),(0.001,0.001),(pi,0)]
    pathticks = ["(π/2,π/2)","(π,0)","(π,π)","(π/2,π/2)","(0,0)","(π,0)"]
    Nk = 75  #75
    k_vec,kticks_positioins = create_brillouin_zone_path(path, Nk)
    =#
    #=
    for i in 1:20
        println(getSitePosition(hte_lattice1.lattice,i+1).-getSitePosition(hte_lattice1.lattice,i),getSitePosition(hte_lattice2.lattice,i+1).-getSitePosition(hte_lattice2.lattice,i))
    end
    =#
    #=
    for k_pos in eachindex(k_vec)
        k = k_vec[k_pos]
        println(get_c_k(k,c_iipDyn_mat1,hte_lattice1).-get_c_k(k,c_iipDyn_mat2,hte_lattice2))
    end
    =#
    #isequal(c_iipDyn_mat1, c_iipDyn_mat2)
end