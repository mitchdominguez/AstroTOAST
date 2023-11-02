using Graphs, MetaGraphs
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                              FRAME CONVERSION FRAMEWORK
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
"""
    frameconvert(state, f1::ReferenceFrame, f2::ReferenceFrame)

Function to convert the `state` from coordinates expressed in frame `f1` to 
coordinates expressed in frame `f2`
"""
function frameconvert(state, f1::ReferenceFrame, f2::ReferenceFrame)

    # Check that neither f1 nor f2 are relative frames
    if isrelativeframe(f1) || isrelativeframe(f2)
        ErrorException("Cannot use only one state to convert into a relative frame") |> throw
    end

    numframes = nv(fc_graph)
    dist_g = gdistances(fc_graph, fc_graph[f1, :frame])[fc_graph[f2, :frame]]
    if dist_g >= numframes
        ErrorException("No such frame transformation exists in AstroTOAST") |> throw
    end
    path = a_star(fc_graph, fc_graph[f1, :frame], fc_graph[f2, :frame])
    dist = length(path)

    if dist == 0
        # f1 and f2 are the same frame, so the state should just be returned as normal
        return state
    elseif dist == 1
        # There exists a function to directly convert between f1 and f2
        return fc(state, f1, f2)
    else
        # Traverse the graph of frame conversions one at a time until arriving
        # at the final state
        newstate = state
        for p in path
            newstate = fc(newstate, fc_graph[src(p), :frame], fc_graph[dst(p), :frame])
        end
        return newstate
    end
end

"""
    frameconvert(target, chaser, f1::ReferenceFrame, f2::ReferenceFrame)

Function to convert the `target` state from coordinates expressed in frame `f1` to 
coordinates expressed in frame `f2`
"""
function frameconvert(target, chaser, f1::ReferenceFrame, f2::ReferenceFrame)
    numframes = nv(fc_graph)
    dist_g = gdistances(fc_graph, fc_graph[f1, :frame])[fc_graph[f2, :frame]]
    if dist_g >= numframes
        ErrorException("No such frame transformation exists in AstroTOAST") |> throw
    end
    path = a_star(fc_graph, fc_graph[f1, :frame], fc_graph[f2, :frame])
    dist = length(path)

    # Check if either f1 or f2 is a relative frame
    f1_isrelative = isrelativeframe(f1)
    f2_isrelative = isrelativeframe(f2)

    if dist == 0 && f1 == f2
        # f1 and f2 are the same frame, so the state should just be returned as normal
        return (target, chaser)
    elseif dist == 1
        # There exists a function to directly convert between f1 and f2
        # return fc(state, f1, f2)

        if !(f1_isrelative) && !(f2_isrelative)
            # Case where both f1 and f2 are not relative
            return (fc(target, f1, f2), fc(chaser, f1, f2))
        else
            # Case where one of f1 or f2 are a relative frame
            return fc(target, chaser, f1, f2)
        end

    else
        # Traverse the graph of frame conversions one at a time until arriving
        # at the final state
        #
        if !(f1_isrelative) && !(f2_isrelative)
            # Case where both f1 and f2 are not relative
            return (frameconvert(target, f1, f2), frameconvert(chaser, f1, f2))
        else
            # Case where at least one of f1 or f2 are a relative frame
            newtarget = target
            newchaser = chaser
            for p in path
                src_f = fc_graph[src(p), :frame]
                dst_f = fc_graph[dst(p), :frame]
                if !(isrelativeframe(src_f)) && !(isrelativeframe(dst_f))
                    # Case where both src_f and dst_f are not relative
                    newtarget, newchaser = (frameconvert(newtarget, src_f, dst_f), frameconvert(newchaser, src_f, dst_f))
                else
                    # Case where at least one of src_f and dst_f are a relative frame
                    newtarget, newchaser = fc(newtarget, newchaser, fc_graph[src(p), :frame], fc_graph[dst(p), :frame])
                end
            end
            return newtarget, newchaser
        end

    end
end

"""
    frameconvert(state, epoch::T, f1::ReferenceFrame, f2::ReferenceFrame) where {T<:Real}

Function to convert the `state` from coordinates expressed in frame `f1` to 
coordinates expressed in frame `f2` at the provided `epoch`
"""
function frameconvert(state, epoch::T, f1::ReferenceFrame, f2::ReferenceFrame) where {T<:Real}
    # Check that neither f1 nor f2 are relative frames
    if isrelativeframe(f1) || isrelativeframe(f2)
        ErrorException("Cannot use only one state to convert into a relative frame") |> throw
    end

    # Call the standard frameconvert method if both f1 and f2 are noninertial
    if !isinertialframe(f1) && !isinertialframe(f2)
        return frameconvert(state, f1, f2)
    end

    numframes = nv(fc_graph)
    dist_g = gdistances(fc_graph, fc_graph[f1, :frame])[fc_graph[f2, :frame]]
    if dist_g >= numframes
        ErrorException("No such frame transformation exists in AstroTOAST") |> throw
    end
    path = a_star(fc_graph, fc_graph[f1, :frame], fc_graph[f2, :frame])
    dist = length(path)

    if dist == 0
        # f1 and f2 are the same frame, so the state should just be returned as normal
        return state
    elseif dist == 1
        # There exists a function to directly convert between f1 and f2
        return fc(state, epoch, f1, f2)
    else
        # Traverse the graph of frame conversions one at a time until arriving
        # at the final state
        newstate = state
        for p in path
            # newstate = fc(newstate, epoch, fc_graph[src(p), :frame], fc_graph[dst(p), :frame])
            newstate = frameconvert(newstate, epoch, fc_graph[src(p), :frame], fc_graph[dst(p), :frame])
        end
        return newstate
    end
    
end


"""
    frameconvert(target, chaser, epoch, f1::ReferenceFrame, f2::ReferenceFrame)

Function to convert the `target` and `chaser` state from coordinates expressed in frame `f1` to 
coordinates expressed in frame `f2` at the provided `epoch`
"""
function frameconvert(target, chaser, epoch::T, f1::ReferenceFrame, f2::ReferenceFrame) where {T<:Real}

    # Call the autonomous frameconvert method if both f1 and f2 are noninertial
    if !isinertialframe(f1) && !isinertialframe(f2)
        return frameconvert(target, chaser, f1, f2)
    end

    numframes = nv(fc_graph)
    dist_g = gdistances(fc_graph, fc_graph[f1, :frame])[fc_graph[f2, :frame]]
    if dist_g >= numframes
        ErrorException("No such frame transformation exists in AstroTOAST") |> throw
    end
    path = a_star(fc_graph, fc_graph[f1, :frame], fc_graph[f2, :frame])
    dist = length(path)

    # Check if either f1 or f2 is a relative frame
    f1_isrelative = isrelativeframe(f1)
    f2_isrelative = isrelativeframe(f2)

    if dist == 0 && f1 == f2
        # f1 and f2 are the same frame, so the state should just be returned as normal
        return (target, chaser)
    elseif dist == 1
        # There exists a function to directly convert between f1 and f2
        # return fc(state, f1, f2)

        if !(f1_isrelative) && !(f2_isrelative)
            # Case where both f1 and f2 are not relative
            return (frameconvert(target, epoch, f1, f2), frameconvert(chaser, epoch, f1, f2))
        else
            # Case where one of f1 or f2 are a relative frame
            return fc(target, chaser, epoch, f1, f2)
        end

    else
        # Traverse the graph of frame conversions one at a time until arriving
        # at the final state
        #
        if !(f1_isrelative) && !(f2_isrelative)
            # Case where both f1 and f2 are not relative
            return (frameconvert(target, epoch, f1, f2), frameconvert(chaser, epoch, f1, f2))
        else

            # Case where at least one of f1 or f2 are a relative frame
            newtarget = target
            newchaser = chaser
            for p in path
                src_f = fc_graph[src(p), :frame]
                dst_f = fc_graph[dst(p), :frame]
                if !(isrelativeframe(src_f)) && !(isrelativeframe(dst_f))
                    # Case where both src_f and dst_f are not relative
                    newtarget, newchaser = (frameconvert(newtarget, epoch, src_f, dst_f), frameconvert(newchaser, epoch, src_f, dst_f))
                else
                    # Case where at least one of src_f and dst_f are a relative frame
                    newtarget, newchaser = fc(newtarget, newchaser, epoch, fc_graph[src(p), :frame], fc_graph[dst(p), :frame])
                end
            end
            return newtarget, newchaser
        end

    end
end



function frameconvert(states::Vector{Vector{T}}, f1::ReferenceFrame, f2::ReferenceFrame) where {T<:Any}
    outstate = similar(states)
    for i = 1:length(states)
        outstate[i] = frameconvert(states[i], f1, f2)
    end
    return outstate
end

function frameconvert(targets::Vector{Vector{T}}, chasers::Vector{Vector{T}}, f1::ReferenceFrame, f2::ReferenceFrame) where {T<:Any}
    if !(length(targets) == length(chasers))
        DimensionMismatch("targets and chasers must be the same length") |> throw
    end

    out_t = similar(targets)
    out_c = similar(chasers)

    for i = 1:length(targets)
        out_t[i], out_c[i] = frameconvert(targets[i], chasers[i], f1, f2)
    end
    return out_t, out_c
end

function frameconvert(states::Vector{Vector{T}}, epochs::AbstractVector, f1::ReferenceFrame, f2::ReferenceFrame) where {T<:Any}
    outstate = similar(states)
    for i = 1:length(states)
        outstate[i] = frameconvert(states[i], epochs[i], f1, f2)
    end
    return outstate
end
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                             REGISTER FRAME CONVERSIONS
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #


# Earth-Moon CR3BP model
em_cr3bp = Cr3bpModel(Bodies["Earth"],Bodies["Moon"])

# Graph to hold frame conversion information
fc_graph = MetaDiGraph()

# Each vertex of fc_graph represents a frame
add_vertex!(fc_graph, :frame, EM_BCR())
add_vertex!(fc_graph, :frame, EM_ECR())
add_vertex!(fc_graph, :frame, EM_MCR())
add_vertex!(fc_graph, :frame, EM_ECAI())
add_vertex!(fc_graph, :frame, EM_MCAI())

# Graph to hold frame conversion information for relative frames
add_vertex!(fc_graph, :frame, EM_TCR())
add_vertex!(fc_graph, :frame, EM_LVLH())
add_vertex!(fc_graph, :frame, EM_VNC())
add_vertex!(fc_graph, :frame, EM_VCR())
add_vertex!(fc_graph, :frame, EM_TCAI())

# Allow each vertex of fc_graph to be indexed by the frame corresponding to it
set_indexing_prop!(fc_graph, :frame)

# Each edge of fc_graph represents an available frame conversion
add_edge!(fc_graph, fc_graph[EM_BCR(), :frame], fc_graph[EM_MCR(), :frame])
add_edge!(fc_graph, fc_graph[EM_MCR(), :frame], fc_graph[EM_BCR(), :frame])

add_edge!(fc_graph, fc_graph[EM_BCR(), :frame], fc_graph[EM_ECR(), :frame])
add_edge!(fc_graph, fc_graph[EM_ECR(), :frame], fc_graph[EM_BCR(), :frame])

add_edge!(fc_graph, fc_graph[EM_BCR(), :frame], fc_graph[EM_TCR(), :frame])
add_edge!(fc_graph, fc_graph[EM_TCR(), :frame], fc_graph[EM_BCR(), :frame])

add_edge!(fc_graph, fc_graph[EM_BCR(), :frame], fc_graph[EM_LVLH(), :frame])
add_edge!(fc_graph, fc_graph[EM_LVLH(), :frame], fc_graph[EM_BCR(), :frame])

add_edge!(fc_graph, fc_graph[EM_BCR(), :frame], fc_graph[EM_VCR(), :frame])
add_edge!(fc_graph, fc_graph[EM_VCR(), :frame], fc_graph[EM_BCR(), :frame])
# add_edge!(fc_graph, fc_graph[EM_BCR(), :frame], fc_graph[EM_VNC(), :frame])

## Arbitrary inertial to the respective frames
add_edge!(fc_graph, fc_graph[EM_MCR(), :frame], fc_graph[EM_MCAI(), :frame])
add_edge!(fc_graph, fc_graph[EM_MCAI(), :frame], fc_graph[EM_MCR(), :frame])
add_edge!(fc_graph, fc_graph[EM_ECR(), :frame], fc_graph[EM_ECAI(), :frame])
add_edge!(fc_graph, fc_graph[EM_ECAI(), :frame], fc_graph[EM_ECR(), :frame])

## Arbitrary inertial absolute to arbitrary inertial target centered
add_edge!(fc_graph, fc_graph[EM_TCAI(), :frame], fc_graph[EM_MCAI(), :frame])
add_edge!(fc_graph, fc_graph[EM_MCAI(), :frame], fc_graph[EM_TCAI(), :frame])
add_edge!(fc_graph, fc_graph[EM_TCAI(), :frame], fc_graph[EM_ECAI(), :frame])
add_edge!(fc_graph, fc_graph[EM_ECAI(), :frame], fc_graph[EM_TCAI(), :frame])
# TODO function to plot the frame conversion graph

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

# """
# `EM_VCR` -> `EM_VNC`
# """
# function fc(target, chaser, f1::EM_VCR, f2::EM_VNC)
    # # SOMETHING ABOUT FLIPPING AROUND THE UNIT VECTORS
# end


### EM_BCR <-> EM_LVLH ###
include("lvlh_frameconversion.jl")

### EM_BCR <-> EM_TCICR ###
include("vcr_frameconversion.jl")

### EM_XCR <-> EM_XCAI
include("arbitrary_inertial_conversion.jl")
