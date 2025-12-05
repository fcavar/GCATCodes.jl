using CairoMakie
using GraphMakie
using Graphs


# -------------------------------------------------- FUNCTIONS --------------------------------------------------


# construct graph from scratch
function construct_graph!(data::CodonGraphData)
    # check if any duplicates in codons
    if length(data.codons) == length(unique(data.codons)) # no duplicates if true
        @debug "node_labels: $(data.node_labels)"
        @debug "length of data.node_labels: $((length(data.node_labels)))"
        @debug "nv(graph): $(nv(data.graph))"
        @debug "nv(graph) == length(data.node_labels): $(nv(data.graph) == length(data.node_labels))"
        @debug "single_bases: $(data.singular_base_nodes)"
        @debug "tuple_base_nodes: $(data.tuple_base_nodes)"
        @debug "node_index: $(data.node_index)"

        create_all_nodes!(data.codons, data.singular_base_nodes, data.tuple_base_nodes)
        create_graph!(data.graph, data.singular_base_nodes, data.tuple_base_nodes, data.node_labels)
        # add an index (id) for each node label
        data.node_index = Dict(label => index for (index, label) in enumerate(data.node_labels))
        @debug "------------------------------"
        @debug "node_labels: $(data.node_labels)"
        @debug "length of data.node_labels: $((length(data.node_labels)))"
        @debug "nv(graph): $(nv(data.graph))"
        @debug "nv(graph) == length(data.node_labels): $(nv(data.graph) == length(data.node_labels))"
        @debug "single_bases: $(data.singular_base_nodes)"
        @debug "tuple_base_nodes: $(data.tuple_base_nodes)"
        @debug "node_index: $(data.node_index)"
        connect_edges!(data.graph, data.codons, data.singular_base_nodes, data.tuple_base_nodes, data.node_index)

        @debug "------------------------------"
        @debug "node_labels: $(data.node_labels)"
        @debug "length of data.node_labels: $((length(data.node_labels)))"
        @debug "nv(graph): $(nv(data.graph))"
        @debug "nv(graph) == length(data.node_labels): $(nv(data.graph) == length(data.node_labels))"
        @debug "single_bases: $(data.singular_base_nodes)"
        @debug "tuple_base_nodes: $(data.tuple_base_nodes)"
        @debug "node_index: $(data.node_index)"
    end
    @debug "Graph constructed from codon set: $(data.codons)"
    show_graph(data)
end


# create all needed nodes for the graph by iterating through codons and collect all singular bases and tuple bases
function create_all_nodes!(codons::Vector{String}, singles::Vector{String}, tuples::Vector{String})
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
function create_graph!(graph::Graphs.DiGraph, singles::Array{String}, tuples::Array{String}, node_lables::Array{String})
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
function connect_edges!(graph::Graphs.DiGraph, codons::Array{String}, singles::Array{String}, tuples::Array{String},
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