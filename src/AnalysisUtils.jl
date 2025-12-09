# -------------------------------------------------- VARIABLES --------------------------------------------------

# -------------------------------------------------- CONSTANTS --------------------------------------------------

# -------------------------------------------------- FUNCTIONS --------------------------------------------------
"""
is_circular(data::CodonGraphData) -> Bool

Returns true if the codon graph represented by `data` is circular aka. is acyclic (does not contain any cycles)
"""
# check if a set of codons is circular by checking if the graph is acyclic
function is_circular(data::CodonGraphData)
    # state_array to keep track of all vertices (0 = unvisited, 1 = visiting, 2 = visited)
    state_array = fill(0, nv(data.graph))

    # perform DFS for each vertice
    for vertice in vertices(data.graph)
        if state_array[vertice] == 0
            if dfs_cycle_detection(data, vertice, state_array)
                return false # no cycles detected, Graph is acyclic aka. circular
            end
        end
    end

    return true # no cycles detected, Graph is not acyclic aka. not circular
end

# recursive DFS to detect cycles
function dfs_cycle_detection(data::CodonGraphData, vertice::Int, state_array::Vector{Int})
    # mark the current vertice as visiting
    state_array[vertice] = 1

    for neighbor in outneighbors(data.graph, vertice)
        if state_array[neighbor] == 1
            return true # neighbor already visited, cycle detected
        elseif state_array[neighbor] == 0
            if dfs_cycle_detection(data, neighbor, state_array)
                return true # cycle detected in recursion
            end
        end
    end

    # mark the current vertice as visited
    state_array[vertice] = 2
    return false # no cycle detected from this path
end


# check if a set of codons is comma-free by checking if a path longer than 2 exists
function is_comma_free(data::CodonGraphData)
    for vertice in vertices(data.graph)
        if dfs_depth_limited(data.graph, vertice, 0)
            return false # path longer than 2 found
        end
    end
    return true # no paths longer than 2 found
end


# recursive depth-limited DFS to find paths longer than 2
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


# check if a set of codons is self-complementary by checking if the graph G(X) is equal to its inverted Graph
function is_self_complementary(data::CodonGraphData)
    # create complement reversed codons
    complement_reversed_codon_set = create_complement_reversed_codons(data)
    # create new graph with complement_reversed_codon_set to compare with original graph
    data_self_complementary = CodonGraphData(
        Graphs.SimpleDiGraph(0),
        complement_reversed_codon_set,
        Vector{String}(),
        Vector{Tuple{String,String}}(),
        Dict{String,Int}(),
    )
    construct_graph!(data_self_complementary)
    return is_graphs_identical(data, data_self_complementary)
end


# compare two graphs for identical structure
function is_graphs_identical(data_first::CodonGraphData, data_second::CodonGraphData)
    # check if same amount of vertices and edges
    @debug "Comparing graphs..."
    @debug "Graph 1: nv=$(nv(data_first.graph)), ne=$(ne(data_first.graph))"
    @debug "Graph 2: nv=$(nv(data_second.graph)), ne=$(ne(data_second.graph))"
    if nv(data_first.graph) != nv(data_second.graph)
        @debug "Not the same amount of vertices"
        return false
    end
    if ne(data_first.graph) != ne(data_second.graph)
        @debug "Not the same amount of edges"
        return false
    end

    # check if same vertice labels
    @debug "Comparing vertice labels..."
    for index in 1:nv(data_first.graph)
        @debug """
        In Graph 1: vertice $(index): $(data_first.vertice_labels[index])
               In Graph 2: vertice $(index): $(data_second.vertice_labels[index])
        """
        if !has_vertice_label(data_second, data_first.vertice_labels[index])
            @debug "vertice label $(data_first.vertice_labels[index]) NOT found in Graph 2"
            return false
        end
    end

    # check if same edges
    @debug "Comparing edges..."
    for edge in edges(data_first.graph)
        src_label = data_first.vertice_labels[src(edge)]
        dst_label = data_first.vertice_labels[dst(edge)]
        @debug "Edge: $(data_first.vertice_labels[src(edge)]) -> $(data_first.vertice_labels[dst(edge)])"
        if has_edge_label(data_second, src_label, dst_label)
            @debug "Edge: $(data_first.vertice_labels[src(edge)]) -> $(data_first.vertice_labels[dst(edge)]) also in Graph 2"
        else
            @debug "Edge: $(data_first.vertice_labels[src(edge)]) -> $(data_first.vertice_labels[dst(edge)]) NOT in Graph 2-------------------------------------"
            return false
        end
    end

    return "Graphs are identical -> self-complementary is true"
    # return "Graph for codon set $(data_first.codon_set) and $(data_second.codon_set) are identical"
end


# check if edge exists in graph between two vertice labels
function has_edge_label(data::CodonGraphData, src_label::String, dst_label::String)
    # check if labels exist
    haskey(data.vertice_index, src_label) || return false
    haskey(data.vertice_index, dst_label) || return false

    start_id = data.vertice_index[src_label]
    end_id = data.vertice_index[dst_label]

    return has_edge(data.graph, start_id, end_id)
end


# check if a vertice labels exists in the graph
function has_vertice_label(data::CodonGraphData, label::String)
    return haskey(data.vertice_index, label)
end