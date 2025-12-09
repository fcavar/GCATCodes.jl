module GCATCodes

# load other source files
include("Types.jl")
include("CodonUtils.jl") # needs Types.jl
include("AnalysisUtils.jl") # needs Types.jl
include("PlotUtils.jl") # needs Types.jl
include("GraphUtils.jl") # needs Types.jl and PlotUtils.jl

# export relevant types
export
# Types.jl
    CodonGraphData
#

# export relevant functions
export
    # AnalysisUtils.jl
    is_circular,
    is_comma_free,
    is_graphs_identical,
    is_self_complementary,
    # CodonUtils.jl
    create_complement_reversed_codons,
    # GCATCodes.jl
    # GraphUtils.jl
    construct_graph!,
    # PlotUtils.jl
    show_graph

end
