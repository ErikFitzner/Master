using JLD2, DelimitedFiles, SimpleWeightedGraphs, Plots, Symbolics
using Profile, ProfileView
include("plotConventions.jl")
include("LatticeGraphs.jl")
include("Embedding.jl")
include("ConvenienceFunctions.jl") 

# TODO add an if statement to distinguish case using one bond type from case with different types --> unique_gG_vec vs GraphG

### load graph evaluations
spin_length = 1/2
n_max = 12

### prepare lattice
lattice_type = "chain"
#Γ,K,M = (0,0), (2*π/3,2*π/sqrt(3)), (0,2*π/sqrt(3))

# Are there J1, J2, J3, J4 interactions?
j1 = true
j2 = true
j3 = false
j4 = false

L = 12

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
        c_iipDyn_mat = get_c_iipDyn_mat(hte_lattice,hte_graphs,C_Dict_vec);
        save_object(fileName_c,c_iipDyn_mat)
    end
end

#test = load_object("CaseStudy/Triangular_Lattice/Triangular_Lattice_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(12) * "_L" * string(12) * ".jld2")
#println("test = $(test[1,1])")

#result = get_TGiip_Matsubara_xpoly(c_iipDyn_mat,1,2,0)
#println("result = $(result)")

#c_iipEqualTime_mat = get_c_iipEqualTime_mat(c_iipDyn_mat)
#result_suscept = sum([get_TGiip_Matsubara_xpoly(c_iipEqualTime_mat, j, 1, 0) for j in 1:length(c_iipDyn_mat)])
#println("uniform susceptebility = $(result_suscept)")

#########################################################################################
###### Dynamic structure factor (DSF) ###################################################
#########################################################################################
if false
a = 1.0
#substitute specific values for x2,x3,x4 to reduce to the previous form of c_iipDyn_mat
fileName_c = "CaseStudy/$(lattice_type)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_a_$(a)_b_$(b)_c_$(c).jld2"
if isfile(fileName_c)
    println("loading "*fileName_c)
    c_iipDyn_mat_subst = load_object(fileName_c)
else
    println("substituting c_iipDyn_mat with a=$(a), b=$(b), c=$(c)")
    c_iipDyn_mat_subst = get_c_iipDyn_mat_subst(c_iipDyn_mat,hte_lattice,a,b,c);
    save_object(fileName_c,c_iipDyn_mat_subst)
end

#### moments and δs of the continued fraction expansion
#(we will calculate the first four moments at fixed k)
k = 0.6*pi            #define fixed k 
x_vec = 0:0.01:4.5   #define temperature range of interest
#Fourier transform the correlation functions at k
c_kDyn_mat = get_c_k([(k,0.0)],c_iipDyn_mat_subst,hte_lattice)[1]
#calculate the moments 
m_vec = get_moments_from_c_kDyn(c_kDyn_mat)
#rescale moments 
m_vec_times_x = [m_vec[i]*Polynomial([0,1]) for i=1:length(m_vec)]

##EXTRAPLOATION OF MOMENTS
#basic pade 
m_vec_extrapolated_pade = []
for m_idx=1:length(m_vec)-2
    push!(m_vec_extrapolated_pade   , extrapolate_series(m_vec[m_idx],"pade",(7-m_idx,7-m_idx)))
end

#pade in u = tanh(f*x) (2 different versions)
f = 0.35   #define the f value (f=0.48 is very fine tuned to give good results)
m_vec_times_x_extrapolated_u_pade1 = []
m_vec_times_x_extrapolated_u_pade2 = []
for m_idx=1:length(m_vec)-2
    push!(m_vec_times_x_extrapolated_u_pade1, extrapolate_series(m_vec_times_x[m_idx],"u_pade",(8-m_idx,7-m_idx,f)))
    push!(m_vec_times_x_extrapolated_u_pade2, extrapolate_series(m_vec_times_x[m_idx],"u_pade",(7-m_idx,8-m_idx,f)))
end

#plot the moments 
plt_m = Plots.plot([0],[0],label="",xlabel=L"x",ylabel=L"x \cdot m_{\mathbf{k},r}(x)/m_{\mathbf{k},r}(0)",title="moments at k="*string(k/pi)*"π",legend=:topleft,xlim=(-0.2,4.5),size=(1.5*aps_width,aps_width))
Plots.plot!(plt_m,[0],[0],label="x bare",color = "grey",linestyle = linestyle_vec[1],linewidth=0.4)
Plots.plot!(plt_m,[0],[0],label="u Padé [7-r,6-r]",color = "grey",linestyle = linestyle_vec[2],alpha =0.5)
Plots.plot!(plt_m,[0],[0],label="u Padé [6-r,7-r]",color = "grey",linestyle = linestyle_vec[3])
#Plots.plot!(plt_m,[0],[0],label="Padé [6-r,6-r]",color = "grey",linestyle = linestyle_vec[4])
for i=1:4
    Plots.plot!(plt_m,[0],[0],label="r="*string(i-1),color = color_vec[i])
    Plots.plot!(plt_m,x_vec[1:180],m_vec_times_x[i].(x_vec[1:180])./(m_vec[i](0)),label = "",alpha= 0.7,color = color_vec[i],linestyle = linestyle_vec[1],linewidth=0.5)
    Plots.plot!(plt_m,x_vec,m_vec_times_x_extrapolated_u_pade1[i].(tanh.(f.*x_vec))./(m_vec[i](0)),label = "",alpha= 0.5,color = color_vec[i],linestyle = linestyle_vec[2],linewidth=0.5)
    Plots.plot!(plt_m,x_vec,m_vec_times_x_extrapolated_u_pade2[i].(tanh.(f.*x_vec))./(m_vec[i](0)),label = "",alpha= 1,color = color_vec[i],linestyle = linestyle_vec[3],linewidth=0.5)
    #Plots.plot!(plt_m,x_vec,x_vec.*m_vec_extrapolated_pade[i].(x_vec)./(m_vec_extrapolated_pade[i](0.0001)),label = "",alpha= 1,color = color_vec[i],linestyle = linestyle_vec[4],linewidth=0.5)
end

#SHOW PLOT
##############
display(plt_m)

#end
#if false

#now calculate the δs of the continued fraction expansion 
#(we will use the u_pade extrapolation)
x_vec = [0.5,1.0,2.0,4.0]  #define the x values at which we calculate the deltas

deltas_at_x =[]
for x in x_vec
    δ_vec,r_vec = fromMomentsToδ([m(tanh(f*x))/x for m in m_vec_times_x_extrapolated_u_pade1])
    push!(deltas_at_x,δ_vec)
end

#plots the δs 
plt_δ = Plots.plot([0],[0],title="δs using Padé with substitution on moments",label="",xlabel=L"r",ylabel=L"\delta_r",legend=:topleft,size=(1.5*aps_width,aps_width))
for (i,x) in enumerate(x_vec)
    Plots.scatter!(plt_δ,0:3,deltas_at_x[i][1:4],label = "x="*string(x),color = thermalCol4_vec[i],markersize=7)
end

#SHOW PLOT
##############
display(plt_δ)




#one can now extrapolate the deltas to infinity via
deltas_at_x_ext = []
for x_idx = 1:length(x_vec)
    δ_vec_ext = extrapolate_δvec(deltas_at_x[x_idx],3,3,1000,false)
    push!(deltas_at_x_ext,δ_vec_ext)
end

#add the extrapolated deltas to the plot 
for (i,x) in enumerate(x_vec)
    Plots.scatter!(plt_δ,4:10,deltas_at_x_ext[i][5:11],label = "",title="δs extrapolated",color = thermalCol4_vec[i],markersize=4,alpha=0.5)
end

#SHOW PLOT
##############
display(plt_δ)




#we can use the deltas to calculate the spin structure factor (still at fixed k)
w_vec = -3:0.01:3  #define w
x_idx = 3          #choose temperature index (relativ to x_vec = [0.5,1.0,2.0,4.0]) 
η = 0.01           #the imaginary part after analytic continuation (η ->0)
#now calculate the DSF 
DSF = [JS(deltas_at_x_ext[x_idx] ,x_vec[x_idx],w,η) for w in w_vec]

plt_dsf = Plots.plot(w_vec,DSF,label="",title="DSF at k="*string(k/pi)*"π and x="*string(x_vec[x_idx]),xlabel=L"\omega/J",ylabel=L"JS(k,\omega)",size = (1.5*aps_width,aps_width))

#SHOW PLOT
##############
display(plt_dsf)

#end
#if true

###SPIN STRUCTURE FACTOR HEATMAPS
using CairoMakie
x = 4.0                #define temperature (x=J/T)
a = 0.5                #define J2/J1
b = 0.0                #define J3/J1
c = 0.0                #define J4/J1

f = 0.6                 # a=0.0->f=0.48, a=0.2->f=0.4, a=0.4->f=0.58, a=0.45->f=0.6, a=0.7-> f=0.65, a=1.0->f=0.35 fits best for x in {0,2,4}
                        # a=0.0->f=0.48, a=0.2->f=0.2, a=0.4->f=0.58, a=0.45->f=0.65, a=0.7-> f=0.25, a=1.0->f=0.4 fits best for x=4

k_step_size = 1/41     #define k step size (in 1/π)
w_step_size = 0.025    #define ω step size (in 1/J)
#define k and ω vectors 
k_vec = vcat(vcat((0.0001,0.0),[(k*pi,0.0) for k in 0:k_step_size:2][2:end-1]),(1.999*pi,0.0)) #[(k,0.0) for k in 0.01:0.0039*2.4:(2*π-0.01)]
w_vec = collect(-3:w_step_size:3.5)

#substitute specific values for x2,x3,x4 to reduce to the previous form of c_iipDyn_mat
fileName_c = "CaseStudy/$(lattice_type)_" * create_spin_string(spin_length) * "_c_iipDyn_nmax" * string(n_max) * "_L" * string(L) * "_a_$(a)_b_$(b)_c_$(c).jld2"
if isfile(fileName_c)
    println("loading "*fileName_c)
    c_iipDyn_mat_subst = load_object(fileName_c)
else
    println("substituting c_iipDyn_mat with a=$(a), b=$(b), c=$(c)")
    c_iipDyn_mat_subst = get_c_iipDyn_mat_subst(c_iipDyn_mat,hte_lattice,a,b,c);
    save_object(fileName_c,c_iipDyn_mat_subst)
end

#println(get_TGiip_Matsubara_xpoly(c_iipDyn_mat_subst, 1, 2, 0)) #test if substitution worked
    
#calculate the spin structure factor for the given k and ω 
JSkw_mat = get_JSkw_mat("u_pade",x,k_vec,w_vec,c_iipDyn_mat_subst,hte_lattice,r_min=3,r_max=3,r_ext=1000,f=f)

#plot the result
fig = Figure(size=(400,400),fontsize=20)
ax=Axis(fig[1,1],limits=(0,2,-3,3),xlabel=L"k/ \pi",ylabel=L"\omega/J_1",title = "x = $x, J₂/J₁ = $a, f = $f",titlesize=20,xlabelsize=20,ylabelsize=20)
hm=CairoMakie.heatmap!(ax,[k[1]/π for k in k_vec],w_vec,JSkw_mat,colormap=:viridis,colorrange=(0.0,0.45),highclip=:white)

#SHOW PLOT
##############
display(fig)
#save("Images/DSF_x_$(x)_J1_J2_$(a)_f_$f.png", fig)

#check sum rule
#println(sum(JSkw_mat)*k_step_size*w_step_size/2)
end

# Substitute specific values
if false
    a = 0.8
    f_num = subvalue(result,a)
    xs1 = range(0, 1.5, length=150)
    xs2 = range(0, 4, length=400)
    ys = [f_num(x) for x in xs1]
    y_pade = robustpade(f_num,6,6).(xs2)
    y_pade2 = robustpade(f_num,5,5).(xs2)
    y_exact = [exact12(x, a) for x in xs2]

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
    Plots.plot!(xs2, y_exact, color="green", alpha=0.7, label="exact")
        
    Plots.plot!(xs2, y_pade_u, color="orange",linestyle=:dash,alpha=0.7,label="u-Padé[6,6] (f=$g)")
    Plots.plot!(xs2, y_pade_u2, color="orange",linestyle=:dashdot,alpha=0.7,label="u-Padé[5,5] (f=$g)")
        
    display(current())
    #savefig("Images/G_12_J2=$(a)J1_u.png")
end