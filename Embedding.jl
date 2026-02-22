#### find embedding factors, adapted to use GraphG type graphs instead of unique_Gg ####

include("Structs.jl")
include("GraphGeneration.jl")
include("vf2_edited.jl") 

""" find if gG is a symmetric graph with respect to switching the two external legs """
function is_symmetric(gG::GraphG)::Bool
    gg = gG.g
    gg_simple = toSimpleGraph(gg)

    edge_relation(b1,b2) = (gg.weights[src(b1),dst(b1)] == gg.weights[src(b2),dst(b2)])

    # finds if there is an isomorphism by only permuting the internal vertices between the graph and the graph with its external vertices flipped.
    #count = count_subgraphisomorph(gg_simple,gg,edge_relation=edge_relation,jL1 = gG.jjp[2],jL2 = gG.jjp[1],jG1 = gG.jjp[1],jG2 = gG.jjp[2])
    count = count_subgraphisomorph(gg,gg,edge_relation=edge_relation,jL1 = gG.jjp[2],jL2 = gG.jjp[1],jG1 = gG.jjp[1],jG2 = gG.jjp[2])
    
    if count >0 
        return true
    else
        return false
    end
end


""" find embedding factor e of GraphG gG in lattice L
- lattice L (needs to be chosen large enough to avoid boundary effects!)
- external site indices jjp=[j,j'] can be [i,i'] (or [i',i] if gG is not symmetric under exchange of i <--> i')
- assumes that the distance j-jp is smaller or equal to the distance of external vertices of gG
"""
function e_fast(LL::SimpleWeightedGraph{Int,Int},j::Int,jp::Int,gG::GraphG, order::Int)::Vector{Int}
    evec = count_subgraphisomorph_by_J2(LL,gG.g,order,jL1 = j,jL2 = jp,jG1 = gG.jjp[1],jG2 = gG.jjp[2])
    return evec
end

""" Calculate the coefficients of (-x)^n for TG_ii'(iν_m) from embedding factors """    
function Calculate_Correlator_fast(L::SimpleWeightedGraph{Int,Int},ext_j1::Int,ext_j2::Int,max_order::Int,hte_graphs::Vector{Vector{GraphG}},C_Dict_vec::Vector{Vector{Vector{Rational{Int128}}}})::Vector{Matrix{Rational{Int128}}}
    # initialize result_array
    result_array = Matrix{Rational{Int128}}[]

    # Convert L to unweighted graph for distance calculation
    L_unweighted = toSimpleGraph(L)
    ext_dist = dijkstra_shortest_paths(L_unweighted, ext_j1).dists[ext_j2]

    #println("($ext_j1,$ext_j2) ext vertices; unweighted distance = $ext_dist")

    for order in 0:max_order
        #println("($ext_j1,$ext_j2) ext vertices; order = $order")
        result_order = zeros(Rational{Int128},order+1,10)
        for (idx, gg) in enumerate(hte_graphs[order+1])
            # Calculate the distance between the external vertices of gg
            gg_unweighted = toSimpleGraph(gg.g)
            gg_dist = dijkstra_shortest_paths(gg_unweighted, gg.jjp[1]).dists[gg.jjp[2]]
            # Apply distance filtering
            if gg_dist < ext_dist
                continue
            end
            if ext_dist == 0
                if gg_dist > ext_dist
                    continue
                end
            else
                if gg_dist == 0
                    continue
                end
            end

            #### Symmetry logic ####
            if is_symmetric(gg) || gg.jjp[1] == gg.jjp[2]
                emb_fac = e_fast(L, ext_j1, ext_j2, gg, order+1)
            else
                emb_fac = e_fast(L, ext_j1, ext_j2, gg, order+1) .+ e_fast(L, ext_j2, ext_j1, gg, order+1)
            end

            #println("($order,$idx) th graph embeding factor = $emb_fac")

            look_up_dict = C_Dict_vec[order+1][idx]
            result_order += emb_fac .* look_up_dict' .* (1//symmetryFactor(gg))
            #println(emb_fac .* look_up_dict' .* (1//symmetryFactor(gg)))
        end
        push!(result_array, result_order)
    end
    return result_array
end