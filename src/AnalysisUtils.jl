# -------------------------------------------------- FUNCTIONS --------------------------------------------------


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