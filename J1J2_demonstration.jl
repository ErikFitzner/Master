using JLD2, DelimitedFiles, SimpleWeightedGraphs, Plots, Symbolics

include("plotConventions.jl")
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 

@variables x1 x2 δ

### prepare lattice ####
spin_length = 1/2
n_max = 12
L = 12
lattice_type = "square"

# TODO properly implement J3,J4 interactions as well, both in "get_finite_Lattice" and in overall logic
# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = true
j3 = false
j4 = false

hte_lattice = getLattice(L,lattice_type,j1,j2,j3,j4);

#### plot lattice: different bond colors for J1,J2,J3,J4 ####
if false
    weights = hte_lattice.graph.weights
    # Assign a color for each bond weight: 1=blue, 2=red, 3=green, 4=orange, else=gray
    color_map = Dict(1.0 => :blue, 2.0 => :red, 3.0 => :green, 4.0 => :orange)
    edge_colors = [get(color_map, w, :gray) for w in weights]
    display(graphplot(hte_lattice.graph, names=1:nv(hte_lattice.graph), markersize=0.2, fontsize=8, nodeshape=:rect, curves=false, edgecolor=edge_colors))
end

#### compute all correlations in the lattice (or load them) ####
if true
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


#########################################################################################
#### Symbolic results: x1=J1/T, x2=J2/T=(J2/J1)*x1 ####
#c_iipDyn_mat_symbolic = get_c_iipDyn_mat_symbolic(c_iipDyn_mat,hte_lattice)

#### test single Matsubara correlation function ####
#i = 1
#iprime = 1
#m = true  # m=0?
#result = get_TGiip_Matsubara_xpoly_symbolic(c_iipDyn_mat_symbolic,i,iprime,m,n_max)
#println("result = $(result)")

#### test uniform susceptibility ####
#c_iipEqualTime_mat = get_c_iipEqualTime_mat_symbolic(c_iipDyn_mat_symbolic)
#result_suscept = sum([get_TGiip_Matsubara_xpoly_symbolic(c_iipEqualTime_mat, j, 1, true, n_max) for j in 1:hte_lattice.lattice.length])
#println("uniform susceptebility = $(expand(result_suscept))")

#### substitute x2 ####
#expr_sub = substitute(result_suscept, Dict(x2 => 0))
#println("uniform susceptibility (x2=0) = $(expr_sub)")
#########################################################################################


#########################################################################################
###### Dynamic structure factor (DSF) ###################################################
#########################################################################################
# Note: For the "chain" geometry, the case a=0 cretes problems 
if true
    #### substitute specific values for x2,x3,x4 to reduce to the previous form of c_iipDyn_mat ####
    a = 1.0 # J_2/J_1
    b = 0.0 # J_3/J_1
    c = 0.0 # J_4/J_1

    sc = sqrt(1+a^2+b^2+c^2) # J/J1 --> scale to w/J and JS
    
    fileName_c = "CaseStudy/$(lattice_type)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_a_$(a)_b_$(b)_c_$(c).jld2"
    if isfile(fileName_c)
        println("loading "*fileName_c)
        c_iipDyn_mat_subst = load_object(fileName_c)
    else
        println("substituting c_iipDyn_mat with a=$(a), b=$(b), c=$(c)")
        c_iipDyn_mat_subst = get_c_iipDyn_mat_subst(c_iipDyn_mat,hte_lattice,a,b,c);
        save_object(fileName_c,c_iipDyn_mat_subst)
    end

    #### check moments ####
    if false
        #(we will calculate the first four moments at fixed k)
        k = pi  #define fixed k 
        x_vec = 0:0.01:2.0  #define temperature range of interest
        #Fourier transform the correlation functions at k
        c_kDyn_mat = get_c_k([(k,k)],c_iipDyn_mat_subst,hte_lattice)[1]
        #calculate the moments 
        m_vec = get_moments_from_c_kDyn(c_kDyn_mat)
        #rescale moments 
        m_vec_times_x = [m_vec[i]*Polynomial([0,1]) for i=1:length(m_vec)]

        ##EXTRAPLOATION OF MOMENTS
        #basic pade 
        m_vec_extrapolated_pade = []
        for m_idx=1:length(m_vec)-2
            push!(m_vec_extrapolated_pade, extrapolate_series(m_vec[m_idx],"pade",(7-m_idx,7-m_idx)))
        end

        #pade in u = tanh(f*x) (2 different versions)
        f = 0.7   #define the f value (f=0.48 is very fine tuned to give good results)
        m_vec_times_x_extrapolated_u_pade1 = []
        m_vec_times_x_extrapolated_u_pade2 = []
        for m_idx=1:length(m_vec)-2
            push!(m_vec_times_x_extrapolated_u_pade1, extrapolate_series(m_vec_times_x[m_idx],"u_pade",(8-m_idx,7-m_idx,f)))
            push!(m_vec_times_x_extrapolated_u_pade2, extrapolate_series(m_vec_times_x[m_idx],"u_pade",(7-m_idx,8-m_idx,f)))
        end

        #plot the moments 
        plt_m = Plots.plot([0],[0],label="",xlabel=L"x",ylabel=L"x \cdot m_{\mathbf{k},r}(x)/m_{\mathbf{k},r}(0)",title="moments at k="*string(k/pi)*"π",legend=:topleft,xlim=(-0.2,2.0),size=(1.5*aps_width,aps_width))
        Plots.plot!(plt_m,[0],[0],label="x bare",color = "grey",linestyle = linestyle_vec[1],linewidth=0.4)
        Plots.plot!(plt_m,[0],[0],label="u Padé [7-r,6-r]",color = "grey",linestyle = linestyle_vec[2],alpha =0.5)
        Plots.plot!(plt_m,[0],[0],label="u Padé [6-r,7-r]",color = "grey",linestyle = linestyle_vec[3])
        #Plots.plot!(plt_m,[0],[0],label="Padé [6-r,6-r]",color = "grey",linestyle = linestyle_vec[4])
        for i=1:4
            Plots.plot!(plt_m,[0],[0],label="r="*string(i-1),color = color_vec[i])
            Plots.plot!(plt_m,x_vec[1:150],m_vec_times_x[i].(x_vec[1:150])./(m_vec[i](0)),label = "",alpha= 0.7,color = color_vec[i],linestyle = linestyle_vec[1],linewidth=0.5)
            Plots.plot!(plt_m,x_vec,m_vec_times_x_extrapolated_u_pade1[i].(tanh.(f.*x_vec))./(m_vec[i](0)),label = "",alpha= 0.5,color = color_vec[i],linestyle = linestyle_vec[2],linewidth=0.5)
            Plots.plot!(plt_m,x_vec,m_vec_times_x_extrapolated_u_pade2[i].(tanh.(f.*x_vec))./(m_vec[i](0)),label = "",alpha= 1,color = color_vec[i],linestyle = linestyle_vec[3],linewidth=0.5)
            #Plots.plot!(plt_m,x_vec,x_vec.*m_vec_extrapolated_pade[i].(x_vec)./(m_vec_extrapolated_pade[i](0.0001)),label = "",alpha= 1,color = color_vec[i],linestyle = linestyle_vec[4],linewidth=0.5)
        end

        #SHOW PLOT
        ##############
        display(plt_m)
    end

    #### compute JS(k,ω) ####
    w_vec = collect(0.0:0.025:3.5*sc)               
    x = 2.0  # J/T
    x0 = x/sc  # J1/T
    f = 0.7
    r_max = 3

    #### define and generate k-path ####
    #path = [(0.001,0.0),(pi,0),(1.999*pi,0)]
    #pathticks = ["0","π","2π"]

    path = [(pi/2,pi/2),(pi,0),(pi,pi),(pi/2,pi/2),(0.001,0.001),(pi,0)]
    pathticks = ["(π/2,π/2)","(π,0)","(π,π)","(π/2,π/2)","(0,0)","(π,0)"]

    #path = [(4*pi/3,0),(0.001,0.001)]
    #pathticks = ["(4π/3,0)","(0,0)"]

    #path = [(0.001,0.001,0.001),(pi,0.001,0.001),(pi,pi,0.001),(pi,pi,pi),(pi/2,pi/2,pi/2),(0.001,0.001,0.001)]
    #pathticks = ["(0,0,0)","(π,0,0)","(π,π,0)","(π,π,π)","(π/2,π/2,π/2)","(0,0,0)"]


    Nk = 75
    k_vec,kticks_positioins = create_brillouin_zone_path(path, Nk)
    
    JSkw_mat = get_JSkw_mat_neu("u_pade",x0,k_vec,w_vec,c_iipDyn_mat_subst,hte_lattice,f=f)
end

#### plot JS(k,ω) ####
if true
    using CairoMakie

    fig = Figure(fontsize=8,size=(aps_width,0.6*aps_width));
    ax=Axis(fig[1,1],xlabel=L"\mathbf{k}",ylabel=L"\omega/J=w",xlabelsize=8,ylabelsize=8);
    hm=CairoMakie.heatmap!(ax,collect(0:Nk)/(Nk),w_vec ./ sc, JSkw_mat .* sc,colormap=:viridis,colorrange=(0.001,0.4),highclip=:white); #0.4
    ax.xticks = ((kticks_positioins .- 1)/(Nk),pathticks)
    CairoMakie.Colorbar(fig[:, end+1], hm,size=8, label = L"J S(\mathbf{k},\omega)")
    CairoMakie.text!(ax,"x=J/T="*string(x),position=[(0.05,0.8)],color=:white)
    CairoMakie.text!(ax,"J2/J1="*string(a),position=[(0.05,0.5)],color=:white)
    CairoMakie.text!(ax,"f=$f",position=[(0.05,0.2)],color=:white)

    resize_to_layout!(fig);
    display(fig)

    #save("Images/$(lattice_type)_JSkw_a_$(a).png",fig; px_per_unit=6.0)
end

#### plot w slice ####
if false
    slice = 81
    wslice = w_vec[slice]
    fig = Figure(fontsize=8,size=(aps_width,0.6*aps_width));
    ax=Axis(fig[1,1],xlabel=L"\mathbf{k}",ylabel=L"JS(k,\omega)",xlabelsize=8,ylabelsize=8, title="$(lattice_type), w=$(wslice)");
    lines!(ax,collect(0:Nk) ./ Nk,JSkw_mat[:,slice].*sc)
    ax.xticks = ((kticks_positioins .- 1) ./ Nk,pathticks)
    display(fig)

    #save("Images/$(lattice_type)_wslice_$(wslice).png",fig; px_per_unit=6.0)
end

#### plot k slice ####
if false
    slice = 40
    kslice = k_vec[slice]
    fig = Figure(fontsize=8,size=(aps_width,0.6*aps_width));
    ax=Axis(fig[1,1],xlabel=L"\omega/J",ylabel=L"JS(k,\omega)",xlabelsize=8,ylabelsize=8, title="$(lattice_type), k=($(round(kslice[1],digits=2)),$(round(kslice[2],digits=2)))");
    lines!(ax,w_vec./sc,JSkw_mat[slice,:].*sc)
    display(fig)

    #save("Images/$(lattice_type)_wslice_$(kslice).png",fig; px_per_unit=6.0)
end
