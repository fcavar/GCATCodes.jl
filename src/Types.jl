using Graphs

export CodonGraphData


mutable struct CodonGraphData
    codons::Vector{String}
    graph::Graphs.SimpleDiGraph
    singular_base_nodes::Vector{String}
    tuple_base_nodes::Vector{String}
    node_labels::Vector{String}
    node_index::Dict{String,Int}
end