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
    Graphs.SimpleDiGraph(0),
    # codon_x0,
    example_codon_set,
    Vector{String}(),
    Vector{Tuple{String, String}}(),
    Dict{String, Int}(),
)
construct_graph!(data; show_plot = true, show_debug = false)
is_circular(data, show_debug = false)
is_comma_free(data, show_debug = false)
is_self_complementary(data, show_plot = false, show_debug = false)
is_c3(data, show_plot = false, show_debug = false)
# "Displaying cycles in graph:", display_cycles(data)
add_vertice_by_label!(data, "C", show_debug = true)
add_vertice_by_label!(data, "TA", show_debug = true)
add_edge_by_label!(data, "T", "AC", show_debug = true)
show_graph(data, show_debug = false)
add_edge_by_label!(data, "AT", "A", show_debug = true)
show_graph(data, show_debug = false)
add_edge_by_label!(data, "A", "TA", show_debug = true)
show_graph(data, show_debug = false)
add_edge_by_label!(data, "TA", "C", show_debug = true)
show_graph(data, show_debug = false)
is_circular(data, show_debug = false)
is_comma_free(data, show_debug = false)
is_self_complementary(data, show_plot = false, show_debug = false)
is_c3(data, show_plot = false, show_debug = false)








# test self complementarity
codon_x0_self_complementary = create_complement_reversed_codons(data::CodonGraphData)
codon_x0_self_complementary_sorted = sort(codon_x0_self_complementary)
println(codon_x0)
println(codon_x0_self_complementary)
println(codon_x0_self_complementary_sorted)
# second data for seconds graph with codon_x0_self_complementary OR codon_x0_self_complementary_sorted
data_self_complementary = CodonGraphData(
    Graphs.SimpleDiGraph(0),
    codon_x0_self_complementary,
    # codon_x0_self_complementary_sorted,
    Vector{String}(),
    Vector{Tuple{String, String}}(),
    Dict{String, Int}(),
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
