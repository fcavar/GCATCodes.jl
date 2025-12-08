using GCATCodes
using CairoMakie
using GraphMakie
using Graphs


# debug logging
using Logging
global_logger(ConsoleLogger(Logging.Debug)) # activate
# global_logger(ConsoleLogger(Logging.Info)) # deactivate


# -------------------------------------------------- FUNCTIONS --------------------------------------------------

# convert a string of codons separated by spaces into a vector of codon strings
function string_to_codon_vector(codon_string::String)
    codon_vector = split(codon_string)
    return codon_vector
end

cod = split("AAC GTT AAG CTT AAT ATT ACC GGT ACG CGT ACT AGT AGC GCT AGG CCT CCG CGG TCA TGA")
@debug "cod: $cod"
@show cod

function test()
    data = CodonGraphData(
        cod,
        Graphs.SimpleDiGraph(0),
        String[],
        String[],
        String[],
        Dict{String,Int}(),
    )
    construct_graph!(data)
end
test()

codon_x0 = ["AAC", "AAT", "ACC", "ATC", "ATT", "CAG", "CTC", "CTG", "GAA", "GAC", "GAG", "GAT", "GCC", "GGC",
    "GGT", "GTA", "GTC", "GTT", "TAC", "TTC"]

data = CodonGraphData(
    codon_x0,
    # ["AAC", "GTT", "AGT", "CGA", "TTC", "GGA", "CTA"],
    Graphs.SimpleDiGraph(0),
    String[],
    String[],
    String[],
    Dict{String,Int}(),
)
construct_graph!(data)
println("Graph is circular: ", is_circular(data))


codon_x0_self_complementary = is_self_complementary(data::CodonGraphData)
data_self_complementary = CodonGraphData(
    codon_x0_self_complementary,
    Graphs.SimpleDiGraph(0),
    String[],
    String[],
    String[],
    Dict{String,Int}(),
)
construct_graph!(data_self_complementary)
println("Graph is circular: ", is_circular(data))


for i in eachindex(codon_x0)
    println("Original codon: $(codon_x0[i]) -> reverse complement codon: $(codon_x0_self_complementary[i])")
end

function is_graphs_identical(g1::Graphs.DiGraph, g2::Graphs.DiGraph)
    nv(g1) == nv(g2) || return false
    ne(g1) == ne(g2) || return false
    for e in edges(g1)
        has_edge(g2, src(e), dst(e)) || return false
    end
    for e in edges(g2)
        has_edge(g1, src(e), dst(e)) || return false
    end

    return true
end

is_graphs_identical(data.graph, Graphs.reverse(data_self_complementary.graph))

isequal(data.graph, Graphs.reverse(data_self_complementary.graph))

# Example plot to test if everything is working
labels = ["A", "B"]
g = Graphs.SimpleDiGraph(2)
Graphs.SimpleGraphs.add_edge!(g, 1, 2)

function a()
    fig = Figure(size=(1800, 900))
    ax = Axis(fig[1, 1])
    ax.title = "Title"
    graphplot!(ax, g;
        nlabels=labels,
        nlabels_color=:white,
        nlabels_size=18,
        nlabels_offset=Point2f(0, 0),
        nlabels_distance_rel=false,
        nlabels_align=(:center, :center),
        node_color=:black,
        node_size=30,
        arrow_shift=:end,
        arrow_size=12,
        edge_width=2,
        edge_curvature=0.9
    )
    fig
end
a()

"""
    demo_function(string::String)

A demo function that takes a string as input and turns it into all capitals.
# Example
```jldoctest
julia> demo_function("hello")
"HELLO"
```
"""

function demo_function(string::String)
    return uppercase(string)
end