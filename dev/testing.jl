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
    codon_x0,
    # example_codon_set,
    Vector{String}(),
    Vector{Tuple{String, String}}(),
    Dict{String, Int}(),
)
construct_graph!(data; show_plot = true, show_debug = false)
println("Graph is circular: ", is_circular(data))
println("Graph is comma-free: ", is_comma_free(data))
println("Graph is self-complementary: ", is_self_complementary(data, show_plot = false))
println("Graph is C3: ", is_c3(data))

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
        construct_graph!(test_data)
        # println(is_self_complementary(test_data))
        if i == 1000
            break
        end
    end
end

v = ["AAC", "AAT"]
function leftshift(s, k = 1)
    k = k % length(s)
    s[k+1:end] * s[1:k]
end
v_shift1 = leftshift.(v)        # um 1 nach links: ["ACA", "ATA"]
v_shift2 = leftshift.(v, 2)     # um 2 nach links: ["CAA", "TAA"]
