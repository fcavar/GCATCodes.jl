using Graphs


# -------------------------------------------------- STRUCTS --------------------------------------------------


# struct to hold all data related to a codon graph
mutable struct CodonGraphData
    graph::Graphs.SimpleDiGraph # directed graph
    codon_set::Vector{String} # the codon set represented in the graph
    vertice_labels::Vector{String} # labels for each vertice
    edge_labels::Vector{Tuple{String,String}} # labels for each edge
    vertice_index::Dict{String,Int} # mapping from vertice label to vertice index in graph
end
