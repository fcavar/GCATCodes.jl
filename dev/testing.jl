isdefined(Main, :GCATCodes) || using GCATCodes
using CairoMakie
using GraphMakie
using Graphs


# debug logging
using Logging
global_logger(ConsoleLogger(Logging.Debug)) # activate
# global_logger(ConsoleLogger(Logging.Info)) # deactivate


# -------------------------------------------------- FUNCTIONS --------------------------------------------------
# create test codon x0 (which is self-complementary, circular and C3)
codon_x0 = ["AAC", "AAT", "ACC", "ATC", "ATT", "CAG", "CTC", "CTG", "GAA", "GAC", "GAG", "GAT", "GCC", "GGC",
    "GGT", "GTA", "GTC", "GTT", "TAC", "TTC"]

# first data for first graph with codon_x0
data = CodonGraphData(
    Graphs.SimpleDiGraph(0),
    codon_x0,
    Vector{String}(),
    Vector{Tuple{String,String}}(),
    Dict{String,Int}(),
)
construct_graph!(data)
println("Graph is circular: ", is_circular(data))

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
    Vector{Tuple{String,String}}(),
    Dict{String,Int}(),
)
construct_graph!(data_self_complementary)
println("Graph is circular: ", is_circular(data))

is_self_complementary(data)

counter = 0
for edge in edges(data.graph)
    counter += 1
    println("Edge $counter in graph 1: $(data.vertice_labels[src(edge)]) -> $(data.vertice_labels[dst(edge)])")
    println("Edge $counter in graph 2: $(data_self_complementary.vertice_labels[src(edge)]) -> $(data_self_complementary.vertice_labels[dst(edge)])")
end


for edge in edges(data.graph)
    if has_edge_labeled(data_self_complementary, data.vertice_labels[src(edge)], data.vertice_labels[dst(edge)])
        println("Edge $(data.vertice_labels[src(edge)]) -> $(data.vertice_labels[dst(edge)]) found in graph 2")
    else
        println("Edge $(data.vertice_labels[src(edge)]) -> $(data.vertice_labels[dst(edge)]) NOT found in graph 2")
    end
end