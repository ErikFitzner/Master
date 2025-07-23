using Graphs, GraphRecipes, Plots
#include Lattice support from SpinMC.jl
include("Structs.jl")
include("Lattice.jl")
include("LatticeSymmetries.jl")

function latticeToGraph(lattice::Lattice)::SimpleWeightedGraph{Int}
    g = SimpleWeightedGraph(lattice.length)
    ints = lattice.interactionSites
    imats = lattice.interactionMatrices

    for v in 1:nv(g)
        for k in 1:length(ints[v])
            u = ints[v][k]
            if u != v && !has_edge(g, v, u)
                M = imats[v][k]
                w = M.m11 ###### Change here for XXZ!
                add_edge!(g, v, u, w)
            end
        end
    end

    return g
end


function get_finite_Lattice(L::Int,geometry::String, j1::Bool, j2::Bool, j3::Bool, j4::Bool; PBC::Bool = true)
    """ creates lattice and corresponding graphs, L is the linear size, PBC sets the use of boundary conditions.
    Currently implmented:
        - chain
        - square
        - triang
        - honeycomb
        - pyrochlore
        - kagome
    """
    # TODO add J2, J3 and J4 interactions everywhere

    a = j1 ? 1 : 0
    b = j2 ? 2 : 0
    c = j3 ? 3 : 0
    d = j4 ? 4 : 0

    M = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    M1 = a*M
    M2 = b*M
    M3 = c*M
    M4 = d*M
    
    if geometry == "chain" ### chain lattice
        a1 = (1, 0)
        a2 = (0, 1)

        uc = UnitCell(a1,a2)
        b0 = addBasisSite!(uc, (0.0,0.0))

        # First neighbors (J1)
        addInteraction!(uc, b0, b0, M1, (1,0))

        # Second neighbors (J2)
        addInteraction!(uc, b0, b0, M2, (2,0))

        # Third neighbors (J3)
        addInteraction!(uc, b0, b0, M3, (3,0))

        # Fourth neighbors (J4)
        addInteraction!(uc, b0, b0, M4, (4,0))

        l = (L,1)

    elseif geometry == "square" ### Square lattice
        a1 = (1, 0)
        a2 = (0, 1)
        uc = UnitCell(a1,a2)

        b0 = addBasisSite!(uc, (0.0, 0.0))

        # First neighbors (J1)
        addInteraction!(uc, b0, b0, M1, (1, 0))
        addInteraction!(uc, b0, b0, M1, (0, 1))

        # Second neighbors (J2, diagonals)
        addInteraction!(uc, b0, b0, M2, (1, 1))
        addInteraction!(uc, b0, b0, M2, (-1, 1))

        # Third neighbors (J3, along axes, distance 2)
        addInteraction!(uc, b0, b0, M3, (2, 0))
        addInteraction!(uc, b0, b0, M3, (0, 2))

        # Fourth neighbors (J4, diagonals at distance sqrt(5))
        addInteraction!(uc, b0, b0, M4, (2, 1))
        addInteraction!(uc, b0, b0, M4, (1, 2))
        addInteraction!(uc, b0, b0, M4, (-2, 1))
        addInteraction!(uc, b0, b0, M4, (-1, 2))

        l = (L, L)

    elseif geometry == "cluster4"
    # Cluster of 4 spins with all-to-all interactions
    a1 = (1, 0)
    a2 = (0, 1)
    uc = UnitCell(a1, a2)

    # Add 4 basis sites
    b0 = addBasisSite!(uc, (0.0, 0.0))
    b1 = addBasisSite!(uc, (1.0, 0.0))
    b2 = addBasisSite!(uc, (0.0, 1.0))
    b3 = addBasisSite!(uc, (1.0, 1.0))

    # Add all-to-all interactions between different spins
    # Non-Diagonals
    addInteraction!(uc, b0, b1, M1, (0, 0))
    addInteraction!(uc, b0, b2, M1, (0, 0))
    addInteraction!(uc, b1, b3, M1, (0, 0))
    addInteraction!(uc, b2, b3, M1, (0, 0))

    # Diagonals
    addInteraction!(uc, b0, b3, M2, (0, 0))
    addInteraction!(uc, b1, b2, M2, (0, 0))

    # Lattice size is just one unit cell — the cluster
    l = (1, 1)


    elseif geometry == "simple_cubic" ### Square lattice
        a1 = (1, 0, 0)
        a2 = (0, 1 ,0)
        a3 = (0, 0, 1)
        uc = UnitCell(a1,a2,a3)

        b0 = addBasisSite!(uc, (0.0, 0.0, 0.0))

        addInteraction!(uc, b0, b0, M1, (1, 0, 0))
        addInteraction!(uc, b0, b0, M1, (0, 1 ,0))
        addInteraction!(uc, b0, b0, M1, (0, 0, 1))

        l = (L, L ,L)

        
    elseif geometry == "triang"  ### Triangular lattice
        a1 = (1/2, sqrt(3)/2)
        a2 = (1/2, -sqrt(3)/2)
        uc = UnitCell(a1,a2)

        b0 = addBasisSite!(uc, (0.0, 0.0))

        addInteraction!(uc, b0, b0, M1, (1, 0))
        addInteraction!(uc, b0, b0, M1, (0, 1))
        addInteraction!(uc, b0, b0, M1, (1, 1))

        #####
        # added for second neighbor interaction, a=J_2/J_1
        addInteraction!(uc, b0, b0, M2, (2, 1))
        addInteraction!(uc, b0, b0, M2, (1, 2))
        addInteraction!(uc, b0, b0, M2, (-1, 1))

        l = (L, L)
 
    elseif geometry == "honeycomb" ### Honeycomb lattice
        a1 = (3/2, sqrt(3)/2)
        a2 = (3/2, -sqrt(3)/2)
        uc = UnitCell(a1,a2)

        b1 = addBasisSite!(uc, (0.0, 0.0))
        b2 = addBasisSite!(uc, (1.0, 0.0))

        addInteraction!(uc, b1, b2, M1, (0, 0))
        addInteraction!(uc, b1, b2, M1, (0, -1))
        addInteraction!(uc, b1, b2, M1, (-1, 0))

        l = (L, L)
        
    elseif geometry == "kagome" ### Honeycomb lattice
        a1 = (1, sqrt(3))
        a2 = (1, -sqrt(3))
        uc = UnitCell(a1,a2)

        b1 = addBasisSite!(uc, (0.0, 0.0))
        b2 = addBasisSite!(uc, -1/2 .* a1)
        b3 = addBasisSite!(uc, 1/2 .* a2)

        addInteraction!(uc, b1, b2, M1, (0, 0))
        addInteraction!(uc, b1, b3, M1, (0, 0))
        addInteraction!(uc, b1, b2, M1, (1, 0))
        addInteraction!(uc, b1, b3, M1, (0, -1))

        addInteraction!(uc, b3, b2, M1, (0, 0))
        addInteraction!(uc, b3, b2, M1, (1, 1))

        

        l = (L, L)

    elseif geometry == "pyrochlore" ### Pyrochlore lattice
        a1 = (0,1/2,1/2)
        a2 = (1/2,0,1/2)
        a3 = (1/2,1/2,0)
        uc = UnitCell(a1,a2,a3)

        b0 = addBasisSite!(uc, (0.0, 0.0, 0.0))
        b1 = addBasisSite!(uc, (0.0, 1/4,1/4))
        b2 = addBasisSite!(uc, (1/4, 0.0 ,1/4))
        b3 = addBasisSite!(uc, (1/4,1/4, 0.0))

        addInteraction!(uc, b0, b1, M1, (0, 0, 0))
        addInteraction!(uc, b0, b2, M1, (0, 0, 0))
        addInteraction!(uc, b0, b3, M1, (0, 0, 0))
        addInteraction!(uc, b0, b1, M1, (-1, 0, 0))
        addInteraction!(uc, b0, b2, M1, (0, -1, 0))
        addInteraction!(uc, b0, b3, M1, (0, 0, -1))

        addInteraction!(uc, b1, b2, M1, (0, 0, 0))
        addInteraction!(uc, b1, b3, M1, (0, 0, 0))
        addInteraction!(uc, b2, b3, M1, (0, 0, 0))
        addInteraction!(uc, b1, b2, M1, (1, -1, 0))
        addInteraction!(uc, b1, b3, M1, (1, 0, -1))
        addInteraction!(uc, b2, b3, M1, (0, 1, -1))

        l = (L, L, L)

    else 
        error("geometry: " * geometry * " not yet implemented") 
    end

    ### create lattice with or without periodic Boundary conditions
    if PBC 
        lattice = Lattice(uc, l);
    else
        lattice = LatticeNoPBC(uc, l);
    end

    ### create graph for the lattice and return
    graph = latticeToGraph(lattice)
    return (lattice,graph)
end

### helper functions
import Base.-
function -(a::Tuple{Float64, Float64},b::Tuple{Float64, Float64})::Tuple{Float64, Float64}
    """ - for 2-tuples """
    return (a[1]-b[1],a[2]-b[2])
end
function -(a::Tuple{Float64, Float64, Float64},b::Tuple{Float64, Float64, Float64})::Tuple{Float64, Float64, Float64}
    """ minus for 3-tuples """
    return (a[1]-b[1],a[2]-b[2],a[3]-b[3])
end

function find_graph_center(graph)
    """ 
    Gives the sites at the center of the graph
    """

    n = nv(graph)  # Number of vertices
    distances = Graphs.floyd_warshall_shortest_paths(graph).dists  # All-pairs shortest path distances
    
    function mean(x)
        sum(x)/length(x)
    end

    # Compute the average distance for each vertex
    avg_distances = [mean(filter(x -> x < typemax(Int), distances[i, :])) for i in 1:n]
    
    # Find the minimum average distance
    min_avg_distance = minimum(avg_distances)
    
    # Find the vertices with the minimum average distance
    min_vertices = findall(x -> x == min_avg_distance, avg_distances)
    
    return min_vertices
end

function getLattice(L::Int,geometry::String, j1::Bool, j2::Bool, j3::Bool,j4::Bool)::Dyn_HTE_Lattice
    """ 
    Gives the lattice where all sites are at most L away from the center sites, no PBC
    """
    # Special case for chain: adjust length for maximum interaction range
    if geometry == "chain"
        # Determine the maximum interaction range present
        max_range = maximum([j1 ? 1 : 0, j2 ? 2 : 0, j3 ? 3 : 0, j4 ? 4 : 0])
        max_range = max(max_range, 1) # Ensure at least 1
        chain_length = 2 * (L * max_range) + 1
        lattice,LatGraph = get_finite_Lattice(chain_length,"chain",j1,j2,j3,j4; PBC = false)
        center_site = div(chain_length + 1, 2)
        center_sites = [center_site]
        lattice.sitePositions = [lattice.sitePositions[i] .- lattice.sitePositions[center_site] for i in 1:chain_length] #shift center site to zero coordinate
        return Dyn_HTE_Lattice(geometry ,lattice, LatGraph, center_sites)
    end
    if geometry == "cluster4"
        lattice, LatGraph = get_finite_Lattice(1, "cluster4",j1,j2,j3,j4; PBC = false)  # assuming you created this case
        center_sites = collect(1:length(lattice.sitePositions))  # all are center
        lattice.sitePositions = [pos .- lattice.sitePositions[1] for pos in lattice.sitePositions]
        return Dyn_HTE_Lattice(geometry ,lattice, LatGraph, center_sites)
    end

    ## helper function
    function replace_indices(tuple_vector,without)
        # Extract all unique integers from the tuples and sort them
        unique_values = unique(collect(Iterators.flatten(tuple_vector)))
        sorted_values = sort([x for x in unique_values if x ∉ without])
        
        # Create a mapping from each value to its ascending order index
        value_to_index = Dict(value => i for (i, value) in enumerate(sorted_values))
        
        # Replace each tuple entry with its corresponding index
        replaced_tuples = [tuple([if x ∉ without value_to_index[x] else i end for x in t ]...) for (i,t) in enumerate(tuple_vector)]
        
        return replaced_tuples
    end


    # Get the lattice and its corresponding graph representation
    lattice, LatGraph = get_finite_Lattice(2 * L + 1,geometry,j1,j2,j3,j4; PBC = false)

    # Extract the number of sites in the unit cell
    basis = length(lattice.unitcell.basis)

    # Compute the indices of the "center" vertices in the lattice
    center_vertices = [basis * sum([(2L + 1)^n * L for n = 0:(length(lattice.size) - 1)]) + b for b in 1:basis]

    # Get the position of the first center vertex (used as a reference point later)
    center_pos = lattice.sitePositions[center_vertices[1]]

    # Calculate the shortest path distances from each center vertex to all other vertices
    # Convert to unweighted graph for distance calculation
    LatGraph_unweighted = toSimpleGraph(LatGraph)
    distances = [dijkstra_shortest_paths(LatGraph_unweighted, center_vertices[i]).dists for i = 1:basis]

    # Identify vertices that are farther away than a threshold distance `L`
    too_large_b = [findall(x -> x > L, distances[i]) for i = 1:basis]

    # Find the common set of "too large" vertices across all basis distances
    too_large = reduce(intersect, [too_large_b[i] for i = 1:basis])

    # Remove the vertices that are too far from the lattice center
    deleteat!(lattice.sitePositions, too_large)
    lattice.length = length(lattice.sitePositions)  # Update the lattice length
    deleteat!(lattice.interactionSites, too_large)  # Remove interactions involving deleted sites

    # Re-center the remaining site positions relative to the initial center vertex
    lattice.sitePositions = [pos - center_pos for pos in lattice.sitePositions]

    # Update the indices of interaction sites after removing "too large" vertices
    lattice.interactionSites = replace_indices(lattice.interactionSites, too_large)

    # Convert the updated lattice back into a graph
    gg = latticeToGraph(lattice)

    # Find the central site(s) of the updated graph (e.g., for further analysis)
    center_sites = find_graph_center(gg)

    # Return the modified lattice, its graph representation, and the central site(s)
    return Dyn_HTE_Lattice(geometry ,lattice, gg, center_sites)
end

###### TESTS ############
#lattice,LatGraph = getLattice(4,"honeycomb"; PBC = false);
#display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))

#lattice,LatGraph,center_sites = getLattice_Ball(6,"honeycomb");
#display(graphplot(LatGraph,names=1:nv(LatGraph),markersize=0.2,fontsize=7,nodeshape=:rect,curves=false))


# Function to determine which basis site a given point corresponds to
function find_site_basis_label(lattice)
    function is_int_vec(vec)
        return all(isinteger,(x->round(x; digits = 10)).(vec))
        end
    
        site_basis_label = zeros(Int,lattice.length)
        basis = collect.(lattice.unitcell.basis)
        lattice_vectors = reduce(hcat,(collect.(lattice.unitcell.primitive)))

        for (site_index,site) in enumerate(collect.(lattice.sitePositions))
        for (i,b) in enumerate(basis)
            fractional_coords = lattice_vectors \ (site - b)
            if is_int_vec(fractional_coords)
                site_basis_label[site_index] = Int(i)
            end
        end
        end

    return site_basis_label
end
