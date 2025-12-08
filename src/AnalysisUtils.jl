# -------------------------------------------------- VARIABLES --------------------------------------------------


# -------------------------------------------------- CONSTANTS --------------------------------------------------
const BASE_COMPLEMENT = Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A')


# -------------------------------------------------- FUNCTIONS --------------------------------------------------

"""
is_circular(data::CodonGraphData) -> Bool

Returns true if the codon graph represented by `data` is circular aka. is acyclic (does not contain any cycles)
"""
# check if a set of codons is circular by checking if the graph is acyclic
function is_circular(data::CodonGraphData)
    # state_array to keep track of all nodes (0 = unvisited, 1 = visiting, 2 = visited)
    state_array = fill(0, nv(data.graph))

    # perform DFS for each node
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
    # mark the current node as visiting
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

    # mark the current node as visited
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
    temp_codons = reverse_codon_set(complement_codon_set(data.codons))
    @debug "Original codon set: $(data.codons) -> reversed complement codon set: $temp_codons"
    return temp_codons
    # inverted_graph = Graphs.reverse(data.graph)
    # return isequal(data.graph, inverted_graph)
end


# returns the complement codon set
function complement_codon_set(codons::Vector{String})
    complement_codons = Vector{String}()
    for codon in codons
        # add the reversed complement codon to the reversed_codons set
        push!(complement_codons, complement_codon(codon))
    end
    @debug "Original codon set: $codons -> complement codon set: $complement_codons"
    return complement_codons
end


# returns the reversed codon set
function reverse_codon_set(codons::Vector{String})
    reversed_codons = Vector{String}()
    for codon in codons
        # add the reversed complement codon to the reversed_codons set
        push!(reversed_codons, reverse_codon(codon))
    end
    @debug "Original codon set: $codons -> reversed codon set: $reversed_codons"
    return reversed_codons
end


# returns the complement base
function complement_base(base::Char)
    @assert base in ['A', 'C', 'G', 'T'] "Base is invalid. Only A, C, G, T are allowed."
    @debug "Original base: $base, -> complement base: $(BASE_COMPLEMENT[base])"
    return BASE_COMPLEMENT[base]
end


# returns the complement codon
function complement_codon(codon::AbstractString)
    @assert all(c in ['A', 'C', 'G', 'T'] for c in codon) "Codon contains invalid characters. Only A, C, G, T are allowed."

    if length(codon) == 3
        @debug "Original codon: $codon, -> complement codon: $(String([complement_base(c) for c in codon]))"
        return String([complement_base(c) for c in codon])
    end
end


# returns the reversed codon
function reverse_codon(codon::AbstractString)
    @assert all(c in ['A', 'C', 'G', 'T'] for c in codon) "Codon contains invalid characters. Only A, C, G, T are allowed."

    if length(codon) == 3
        reversed = reverse(codon)
        @debug "Original codon: $codon, -> reversed codon: $reversed"
        return reversed
    end
end