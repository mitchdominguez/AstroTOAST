using Graphs, MetaGraphs
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                              FRAME CONVERSION FRAMEWORK
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

function frameconvert(state, f1::ReferenceFrame, f2::ReferenceFrame)
    # dist = gdistances(fc_graph, fc_graph[f1, :frame])[fc_graph[f2, :frame]]
    path = a_star(fc_graph, fc_graph[f1, :frame], fc_graph[f2, :frame])
    dist = length(path)
    numframes = length(vertices(fc_graph))

    if dist == 0
        # f1 and f2 are the same frame, so the state should just be returned as normal
        return state
    elseif dist == 1
        # There exists a function to directly convert between f1 and f2
        return fc(state, f1, f2)
    elseif dist < numframes
        # Traverse the graph of frame conversions one at a time until arriving
        # at the final state
        newstate = state
        for p in path
            newstate = fc(newstate, fc_graph[src(p), :frame], fc_graph[dst(p), :frame])
        end
        return newstate
    else
        ErrorException("No such frame transformation exists in AstroTOAST") |> throw
    end
    # return state
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
# add_vertex!(fc_graph, :frame, EM_TCR())
# add_vertex!(fc_graph, :frame, EM_LVLH())
# add_vertex!(fc_graph, :frame, EM_ICR())

# Allow each vertex of fc_graph to be indexed by the frame corresponding to it
set_indexing_prop!(fc_graph, :frame)

# Each edge of fc_graph represents an available frame conversion
add_edge!(fc_graph, fc_graph[EM_BCR(), :frame], fc_graph[EM_MCR(), :frame])
add_edge!(fc_graph, fc_graph[EM_BCR(), :frame], fc_graph[EM_ECR(), :frame])

#TODO TCR

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

"""
`EM_BCR` <-> `EM_TCR`
"""
function fc(target, chaser, f1::EM_BCR, f2::EM_TCR)
    return (target, chaser-target)
end
function fc(target, chaser, f1::EM_TCR, f2::EM_BCR)
    return (target, chaser+target)
end
