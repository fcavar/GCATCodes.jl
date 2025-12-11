using Revise
isdefined(Main, :GCATCodes) || using GCATCodes
using JuliaFormatter
using Logging
using CairoMakie
using GraphMakie
using Graphs


# debug logging
global_logger(ConsoleLogger(Logging.Debug)) # activate
# global_logger(ConsoleLogger(Logging.Info)) # deactivate


# -------------------------------------------------- FUNCTIONS --------------------------------------------------
# create test codon x0 (which is self-complementary, circular and C3)
codon_x0 = [
    "AAC",
    "AAT",
    "ACC",
    "ATC",
    "ATT",
    "CAG",
    "CTC",
    "CTG",
    "GAA",
    "GAC",
    "GAG",
    "GAT",
    "GCC",
    "GGC",
    "GGT",
    "GTA",
    "GTC",
    "GTT",
    "TAC",
    "TTC",
]
example_codon_set = ["CGT", "GTA", "ACT", "AAT"]
# first data for first graph with codon_x0
data = CodonGraphData(
    Graphs.SimpleDiGraph(0), # graph
    codon_x0, # codon_set
    # example_codon_set, # codon_set
    Vector{String}(), # vertice_labels
    Vector{String}(), # added_vertice_labels
    Vector{Tuple{String, String}}(), # edge_labels
    Vector{Tuple{String, String}}(), # added_edge_labels
    Dict{String, Int}(), # vertice_index
)
construct_graph!(data; show_plot = true, show_debug = false)
show_graph(data; show_debug = false)
# check properties of graph
is_circular(data, show_debug = false)
is_comma_free(data, show_debug = false)
is_self_complementary(data, show_plot = false, show_debug = false)
is_c3(data, show_plot = false, show_debug = false)



# manually add vertices and edges to data_adjusted
data_adjusted = CodonGraphData(
    Graphs.SimpleDiGraph(0), # graph
    example_codon_set, # codon_set
    Vector{String}(), # vertice_labels
    Vector{String}(), # added_vertice_labels
    Vector{Tuple{String, String}}(), # edge_labels
    Vector{Tuple{String, String}}(), # added_edge_labels
    Dict{String, Int}(), # vertice_index
)
construct_graph!(data_adjusted; show_plot = true, show_debug = false)

# get get N₂ and N₃N₁ for each codon and add them as vertices and edges between them
println("data_adjusted.codon_set: $(data.codon_set)")
for codon in data_adjusted.codon_set
    n2 = codon[2:2]
    n3n1 = string(codon[3], codon[1])
    println("n2: $n2, n3n1: $n3n1")
    add_vertice_by_label!(data_adjusted, n2, show_debug = false)
    add_vertice_by_label!(data_adjusted, n3n1, show_debug = false)
    connect_edge_by_label!(data_adjusted, n2, n3n1, show_debug = false)
end
show_graph(data_adjusted; show_debug = false)

# check properties of adjusted graph
is_circular(data_adjusted, show_debug = false)
is_comma_free(data_adjusted, show_debug = false)
is_self_complementary(data_adjusted, show_plot = false, show_debug = false)
is_c3(data_adjusted, show_plot = false, show_debug = false)




println(data.codon_set)
println(data_adjusted.codon_set)
println(data.vertice_labels)
println(data_adjusted.vertice_labels)
println("Edges in data")
for edge in edges(data.graph)
    print("$(data.vertice_labels[src(edge)]) -> $(data.vertice_labels[dst(edge)]), ")
end
println("Edges in data.adjusted")
for edge in edges(data_adjusted.graph)
    print(
        "$(data_adjusted.vertice_labels[src(edge)]) -> $(data_adjusted.vertice_labels[dst(edge)]), ",
    )
end











# test self complementarity
codon_x0_self_complementary = create_complement_reversed_codons(data::CodonGraphData)
codon_x0_self_complementary_sorted = sort(codon_x0_self_complementary)
println(codon_x0)
println(codon_x0_self_complementary)
println(codon_x0_self_complementary_sorted)
# second data for seconds graph with codon_x0_self_complementary OR codon_x0_self_complementary_sorted
data_self_complementary = CodonGraphData(
    Graphs.SimpleDiGraph(0), # graph
    codon_x0_self_complementary, # codon_set
    Vector{String}(), # vertice_labels
    Vector{String}(), # added_vertice_labels
    Vector{Tuple{String, String}}(), # edge_labels
    Vector{Tuple{String, String}}(), # added_edge_labels
    Dict{String, Int}(), # vertice_index
)
construct_graph!(data_self_complementary)
println("Graph is circular: ", is_circular(data))

is_self_complementary(data)

counter = 0
for edge in edges(data.graph)
    counter += 1
    println(
        "Edge $counter in graph 1: $(data.vertice_labels[src(edge)]) -> $(data.vertice_labels[dst(edge)])",
    )
    println(
        "Edge $counter in graph 2: $(data_self_complementary.vertice_labels[src(edge)]) -> $(data_self_complementary.vertice_labels[dst(edge)])",
    )
end


for edge in edges(data.graph)
    if has_edge_labeled(
        data_self_complementary,
        data.vertice_labels[src(edge)],
        data.vertice_labels[dst(edge)],
    )
        println(
            "Edge $(data.vertice_labels[src(edge)]) -> $(data.vertice_labels[dst(edge)]) found in graph 2",
        )
    else
        println(
            "Edge $(data.vertice_labels[src(edge)]) -> $(data.vertice_labels[dst(edge)]) NOT found in graph 2",
        )
    end
end


# read file per line with all 216 maximal self-complementary C3 codon_set
open("files/216_maximal_self_complementary_c3_codes.txt", "r") do f
    for (i, line) in enumerate(eachline(f))
        test_data = CodonGraphData(
            Graphs.SimpleDiGraph(0),
            split(line),
            Vector{String}(),
            Vector{Tuple{String, String}}(),
            Dict{String, Int}(),
        )
        println("""Testing codon set #$i:
        $(test_data.codon_set)
        """)
        construct_graph!(test_data, show_plot = false, show_debug = false)
        is_circular(test_data, show_debug = false)
        is_comma_free(test_data, show_debug = false)
        is_self_complementary(test_data, show_plot = false, show_debug = false)
        is_c3(test_data, show_plot = false, show_debug = false)
        if i == 1
            break
        end
    end
end


# example for cycle detection
example_codon_set2 = ["CGT", "GTA", "ACT", "AAT", "ATT", "TTA", "TTC"]
example_data = CodonGraphData(
    Graphs.SimpleDiGraph(0), # graph
    example_codon_set, # codon_set
    Vector{String}(), # vertice_labels
    Vector{String}(), # added_vertice_labels
    Vector{Tuple{String, String}}(), # edge_labels
    Vector{Tuple{String, String}}(), # added_edge_labels
    Dict{String, Int}(), # vertice_index
)
construct_graph!(example_data; show_plot = true, show_debug = false)
show_graph(example_data; show_debug = false)
display_cycles(example_data; show_debug = true)

reverse_data = CodonGraphData(
    Graphs.SimpleDiGraph(0), # graph
    get_reverse_codon_set(example_codon_set), # codon_set
    Vector{String}(), # vertice_labels
    Vector{String}(), # added_vertice_labels
    Vector{Tuple{String, String}}(), # edge_labels
    Vector{Tuple{String, String}}(), # added_edge_labels
    Dict{String, Int}(), # vertice_index
)
construct_graph!(reverse_data; show_plot = true, show_debug = false)

println(example_data.vertice_labels)
println(reverse_data.vertice_labels)

data_list = [example_data, reverse_data, example_data, data, data_adjusted]
show_multiple_graphs(data_list; show_debug = true)