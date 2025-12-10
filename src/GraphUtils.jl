using CairoMakie
using GraphMakie
using Graphs
# -------------------------------------------------- VARIABLES --------------------------------------------------

# -------------------------------------------------- CONSTANTS --------------------------------------------------

# -------------------------------------------------- FUNCTIONS --------------------------------------------------
# construct graph from data
function construct_graph!(data::CodonGraphData; show_plot::Bool = false, show_debug::Bool = false)
    # check if any duplicates in codons
    if length(data.codon_set) == length(unique(data.codon_set)) # no duplicates if true
        show_debug && @debug """
        Before adding vertices and edges:
        graph: $(data.graph)
        codon_set: $(data.codon_set)
        vertice_labels: $(data.vertice_labels)
        edge_labels: $(data.edge_labels)
        """

        # extract vertice labels from codon set and add vertices to graph
        create_all_vertices!(data)
        # create mapping from vertice label to vertice index in graph
        data.vertice_index =
            Dict(label => index for (index, label) in enumerate(data.vertice_labels))
        # connect edges
        connect_edges!(data)

        show_debug && @debug """
        After adding vertices and edges:
        graph: $(data.graph)
        codon_set: $(data.codon_set)
        vertice_labels: $(data.vertice_labels)
        edge_labels: $(data.edge_labels)
        """
    end
    show_debug && @debug "Graph construction from codon set finished: $(data.codon_set)"
    if show_plot
        show_graph(data, show_debug = show_debug)
    end
end


# create all needed vertices for the graph by iterating through codons and collect all singular and tuple bases
function create_all_vertices!(data::CodonGraphData; show_debug::Bool = false)
    # use a temporary set to avoid duplicates and increase lookup speed
    temp_labels = Set{String}()

    for codon in data.codon_set
        # get first and last character of codon
        push!(temp_labels, string(codon[1])) # first base
        push!(temp_labels, string(codon[3])) # third base
        # get first two and last two characters of codon
        push!(temp_labels, codon[1:2]) # first tuple
        push!(temp_labels, codon[2:3]) # second tuple
    end

    # copy Set to vertice_labels
    data.vertice_labels = sort(collect(temp_labels))
    # add vertices to graph
    for _ in 1:length(data.vertice_labels)
        add_vertex!(data.graph)
    end
end


# connect the vertices in the graph based on the codons by connecting the first base to the second tuple and
# the first tuple to the third base
function connect_edges!(data::CodonGraphData; show_debug::Bool = false)
    for codon::AbstractString in data.codon_set
        # get needed vertice IDs
        first_base_id = data.vertice_index[SubString(codon, 1, 1)]
        third_base_id = data.vertice_index[SubString(codon, 3, 3)]
        first_tuple_id = data.vertice_index[SubString(codon, 1, 2)]
        second_tuple_id = data.vertice_index[SubString(codon, 2, 3)]
        add_edge!(data.graph, first_base_id, second_tuple_id)
        add_edge!(data.graph, first_tuple_id, third_base_id)
    end
end
