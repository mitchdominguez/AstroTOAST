using Graphs, MetaGraphs
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                              FRAME CONVERSION FRAMEWORK
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

function frameconvert(state, f1::ReferenceFrame, f2::ReferenceFrame)
    show(stdout, "text/plain", fc_graph)
    return state
end

# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                             REGISTER FRAME CONVERSIONS
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #


# Earth-Moon CR3BP model
em_cr3bp = Cr3bpModel(Bodies["Earth"],Bodies["Moon"])

# Graph to hold frame conversion information
fc_graph = MetaGraph()

# Each vertex of fc_graph represents a frame
add_vertex!(fc_graph, :frame, EM_BCR())
add_vertex!(fc_graph, :frame, EM_ECR())
add_vertex!(fc_graph, :frame, EM_MCR())
add_vertex!(fc_graph, :frame, EM_ECI())
add_vertex!(fc_graph, :frame, EM_MCI())
add_vertex!(fc_graph, :frame, EM_TCR())
add_vertex!(fc_graph, :frame, EM_LVLH())
add_vertex!(fc_graph, :frame, EM_ICR())

# Allow each vertex of fc_graph to be indexed by the frame corresponding to it
set_indexing_prop!(fc_graph, :frame)

# Each edge of fc_graph represents an available frame conversion
add_edge!(fc_graph, fc_graph[EM_BCR(), :frame], fc_graph[EM_MCR(), :frame])
add_edge!(fc_graph, fc_graph[EM_BCR(), :frame], fc_graph[EM_ECR(), :frame])

# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                            FRAME CONVERSION FUNCTIONS
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
`EM_BCR` <-> `EM_MCR`
"""
function fc(state, f1::EM_BCR, f2::EM_MCR)
    return state.-[1-mass_ratio(em_cr3bp), 0, 0, 0, 0, 0]
end
function fc(state, f1::EM_MCR, f2::EM_BCR)
    return state.+[1-mass_ratio(em_cr3bp), 0, 0, 0, 0, 0]
end

"""
`EM_BCR` <-> `EM_ECR`
"""
function fc(state, f1::EM_BCR, f2::EM_ECR)
    return state.+[mass_ratio(em_cr3bp), 0, 0, 0, 0, 0]
end
function fc(state, f1::EM_ECR, f2::EM_BCR)
    return state.-[mass_ratio(em_cr3bp), 0, 0, 0, 0, 0]
end
