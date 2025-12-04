module GCATCodes

# load other source files
include("CodonComplement.jl")
include("Types.jl")
include("PlotUtils.jl") # needs Types.jl
include("GraphUtils.jl") # needs Types.jl and PlotUtils.jl

# write a function that turns a string like "AAC GTT AAG CTT AAT ATT ACC GGT ACG CGT ACT AGT AGC GCT AGG CCT CCG CGG TCA TGA"
# into a vector of codons like ["AAC", "GTT", "AAG",...]
function string_to_codon_vector(codon_string::String)
    codon_vector = split(codon_string)
    return codon_vector
end

cod = split("AAC GTT AAG CTT AAT ATT ACC GGT ACG CGT ACT AGT AGC GCT AGG CCT CCG CGG TCA TGA")
println("cod: $cod")

function test()
    data = CodonGraphData(
        cod,
        Graphs.SimpleDiGraph(0),
        String[],
        String[],
        String[],
        Dict{String,Int}(),
    )
    construct_graph(data)
end
test()


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

end
