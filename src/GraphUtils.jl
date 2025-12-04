using CairoMakie
using GraphMakie
using Graphs

export construct_graph, create_all_nodes, create_graph, show_graph, is_comma_free

# construct graph from scratch
function construct_graph(data::CodonGraphData)
    # check if any duplicates in codons
    if length(data.codons) == length(unique(data.codons)) # no duplicates if true
        println("node_labels: $(data.node_labels)")
        println("length of data.node_labels: $((length(data.node_labels)))")
        println("nv(graph): $(nv(data.graph))")
        println("nv(graph) == length(data.node_labels): $(nv(data.graph) == length(data.node_labels))")
        println("single_bases: $(data.singular_base_nodes)")
        println("tuple_base_nodes: $(data.tuple_base_nodes)")
        println("node_index: $(data.node_index)")

        create_all_nodes(data.codons, data.singular_base_nodes, data.tuple_base_nodes)
        create_graph(data.graph, data.singular_base_nodes, data.tuple_base_nodes, data.node_labels)
        data.node_index = Dict(n => i for (i, n) in enumerate(data.node_labels)) # add an index (id) for each node label        println("------------------------------")
        println("node_labels: $(data.node_labels)")
        println("length of data.node_labels: $((length(data.node_labels)))")
        println("nv(graph): $(nv(data.graph))")
        println("nv(graph) == length(data.node_labels): $(nv(data.graph) == length(data.node_labels))")
        println("single_bases: $(data.singular_base_nodes)")
        println("tuple_base_nodes: $(data.tuple_base_nodes)")
        println("node_index: $(data.node_index)")
        println("---------------")
        println("---------------")
        println("---------------")
        println("---------------")
        println("---------------")
        println("---------------")
        connect_edges(data.graph, data.codons, data.singular_base_nodes, data.tuple_base_nodes, data.node_index)

        println("------------------------------")
        println("node_labels: $(data.node_labels)")
        println("length of data.node_labels: $((length(data.node_labels)))")
        println("nv(graph): $(nv(data.graph))")
        println("nv(graph) == length(data.node_labels): $(nv(data.graph) == length(data.node_labels))")
        println("single_bases: $(data.singular_base_nodes)")
        println("tuple_base_nodes: $(data.tuple_base_nodes)")
        println("node_index: $(data.node_index)")
    end
    println("Graph constructed from codon set: $(data.codons)")
    show_graph(data)
    # show_graph(data.graph, data.node_labels, data.codons)
    # return graph, data.node_labels
end


# create all needed nodes for the graph by iterating through codons and collect all singular bases and tuple bases
function create_all_nodes(codons::Vector{String}, singles::Vector{String}, tuples::Vector{String})
    for codon in codons
        # get first character of codon if not already in "singles"
        if !(string(codon[1]) in singles)
            push!(singles, string(codon[1]))
        end
        # get third character of codon if not already in "singles"
        if !(string(codon[3]) in singles)
            push!(singles, string(codon[3]))
        end
        # get first tuple of codon if not already in "tuples"
        if !(codon[1:2] in tuples)
            push!(tuples, codon[1:2])
        end
        # get second tuple of codon if not already in "tuples"
        if !(codon[2:3] in tuples)
            push!(tuples, codon[2:3])
        end
    end
end


# add all created nodes to the graph and label them
function create_graph(graph::Graphs.DiGraph, singles::Array{String}, tuples::Array{String}, node_lables::Array{String})
    # add all single bases as nodes and label them
    for base::AbstractString in singles
        Graphs.SimpleGraphs.add_vertex!(graph)
        push!(node_lables, base)
    end
    # add all tuple bases as nodes and label them
    for tuple::AbstractString in tuples
        Graphs.SimpleGraphs.add_vertex!(graph)
        push!(node_lables, tuple)
    end

end


# connect the nodes in the graph based on the codons
function connect_edges(graph::Graphs.DiGraph, codons::Array{String}, singles::Array{String}, tuples::Array{String},
    node_index::Dict{String,Int})
    for codon::AbstractString in codons
        # connect first base to second tuple
        first_base_id = node_index[string(codon[1])]
        second_tuple_id = node_index[codon[2:3]]
        Graphs.SimpleGraphs.add_edge!(graph, first_base_id, second_tuple_id)
        # connect first tuple to third base
        first_tuple_id = node_index[codon[1:2]]
        third_base_id = node_index[string(codon[3])]
        Graphs.SimpleGraphs.add_edge!(graph, first_tuple_id, third_base_id)
    end
end


# show graph
function show_graph(graph::Graphs.DiGraph, node_labels::Vector{String}, codons::Vector{String})
    fig = Figure()
    ax = Axis(fig[1, 1])
    ax.title = "Graph for codon set: $codons"
    graphplot!(ax, graph;
        node_labels=node_labels,
        textcolor=:black,
        node_size=30,
        arrow_shift=:end,
        arrow_size=12,
        edge_width=2,
        edge_curvature=0.9
    )
end

# check if a set of codons is comma-free by checking if a path longer than 2 exists
function is_comma_free(graph::Graphs.DiGraph)
    for vertice in vertices(graph)
        if dfs_depth_limited(graph, vertice, 0)
            return false # path longer than 2 found
        end
    end
    return true # no paths longer than 2 found
end


# check if a set of codons is comma-free by checking if a path longer than 2 exists
function dfs_depth_limited(graph::Graphs.DiGraph, vertice::Int, depth::Int)
    if depth >= 3
        return true # Pfad mit Länge >= 3 gefunden
    end

    for neighbor in outneighbors(graph, vertice)
        if dfs_depth_limited(graph, neighbor, depth + 1)
            return true
        end
    end
    return false
end