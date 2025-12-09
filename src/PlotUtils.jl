using CairoMakie
using GraphMakie
# -------------------------------------------------- VARIABLES --------------------------------------------------

# -------------------------------------------------- CONSTANTS --------------------------------------------------

# -------------------------------------------------- FUNCTIONS --------------------------------------------------
# creates a plot representing the codon set data
function show_graph(data::CodonGraphData)
    @debug "Showing graph..."
    fig = Figure(size=(1800, 900))
    ax = Axis(fig[1, 1])
    ax.title = "Graph for codon set: $(data.codon_set)"
    @debug "vertice_labels in graph: $(data.vertice_labels)"
    graphplot!(ax, data.graph;
        nlabels=data.vertice_labels,
        nlabels_color=:white,
        nlabels_size=18,
        nlabels_offset=Point2f(0, 0),
        nlabels_align=(:center, :center),
        node_color=:black,
        node_size=50,
        arrow_shift=:end,
        arrow_size=12,
        edge_width=2,
        edge_curvature=0.9
    )
    fig
end
