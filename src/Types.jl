using Graphs


# -------------------------------------------------- STRUCTS --------------------------------------------------


# struct to hold all data related to a codon graph
mutable struct CodonGraphData
    codons::Vector{String}
    graph::Graphs.SimpleDiGraph
    singular_base_nodes::Vector{String}
    tuple_base_nodes::Vector{String}
    node_labels::Vector{String}
    node_index::Dict{String,Int}
end