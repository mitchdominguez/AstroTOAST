# TODO export initial conditions for each segment with TOFs
using OrdinaryDiffEq



# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                     TRAJECTORY
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
    struct Trajectory{D}

Object that will represent trajectories within a given dynamical model
"""

struct Trajectory{D}
    X_0::Vector{Float64}
    tspan::Vector{Float64}
    dm::DynamicalModel 
    X::Vector{OrdinaryDiffEq.ODESolution} # State history

    # Constructor for generating a trajectory from a single patch point
    function Trajectory(dm::DynamicalModel, X0::AbstractVector{Float64}, proptime)
        # Check that dimensions all match
        if dimension(dm) != length(X0)
            throw(DimensionMismatch("Initial conditions and dynamical model have different dimensions"))
        end

        
        # Create tspan 
        if length(proptime) == 1
            tspan = Vector{Float64}([0, proptime[1]])
        elseif length(proptime) == 2
            tspan = Vector{Float64}([(proptime[1]), proptime[2]])
        else
            throw(DimensionMismatch("proptime must be of length 1 or 2"))
        end

        # Propagate initial conditions for the desired time
        sol = solve(dm, SVector(X0...), tspan)

        new{dimension(dm)}(copy(X0), sort(tspan), dm, [sol])
    end

end

"""
    Trajectory(dm::DynamicalModel, X0, T)

Constructor for generating a trajectory from multiple patch points
"""
function Trajectory(dm::DynamicalModel, X0, T)
    if (length(X0)!=length(T)) && (length(X0)+1!=length(T))
        throw(DimensionMismatch("Different numbers of patch points and times given"))
    end

    lenX = length(X0)
    lenT = length(T)

    if lenX == lenT
        traj = Trajectory(dm, X0[1], T[1])
        for i = 2:lenX
            append!(traj, Trajectory(dm, X0[i], T[i]))
        end
    else
        traj = Trajectory(dm, X0[1], T[1], T[2])
        for i = 2:lenX
            append!(traj, Trajectory(dm, X0[i], T[i], T[i+1]))
        end
    end

    return traj
end

"""
    Trajectory(dm::DynamicalModel, X0::FreeVariable, T::FreeVariable)

Constructor for Trajectory object that uses FreeVariables as inputs
"""
function Trajectory(dm::DynamicalModel, X0::FreeVariable, T::FreeVariable)
    if full_length(X0) == dimension(dm) && full_length(T) == 1
        return Trajectory(dm, tofullvector(X0), tofullvector(T))
    else
        throw(DimensionMismatch("Initial conditions and dynamical model have different dimensions"))
    end
end

"""
    Trajectory(dm::DynamicalModel, X0::FreeVariable, T1::FreeVariable, T2::FreeVariable)

Constructor for Trajectory object that uses FreeVariables as inputs
"""
function Trajectory(dm::DynamicalModel, X0::FreeVariable, T1::FreeVariable, T2::FreeVariable)
    if full_length(X0) == dimension(dm) && full_length(T1) == 1 && full_length(T1) == 1
        return Trajectory(dm, tofullvector(X0), [tofullvector(T1)[1], tofullvector(T2)[1]])
    else
        throw(DimensionMismatch("Initial conditions and dynamical model have different dimensions"))
    end
end

"""
    Trajectory(dm::DynamicalModel, X0::AbstractVector, T1::Float64, T2::Float64)

Constructor for Trajectory object that uses FreeVariables as inputs
"""
function Trajectory(dm::DynamicalModel, X0::AbstractVector, T1::Real, T2::Real)
    if length(X0) == dimension(dm)
        return Trajectory(dm, X0, [T1, T2])
    else
        throw(DimensionMismatch("Initial conditions and dynamical model have different dimensions"))
    end
end

"""
    Trajectory(dm::DynamicalModel, X0::AbstractVector, T1::Float64, T2::Float64)

Constructor for Trajectory object that uses FreeVariables as inputs
"""
function Trajectory(dm::DynamicalModel, X0, T1, T2)
    if (length(X0)!=length(T1)) && length(X0)!=length(T2)
        throw(DimensionMismatch("Different numbers of patch points and times given"))
    end

    lenX = length(X0)

    traj = Trajectory(dm, X0[1], T1[1], T2[1])
    for i = 2:lenX
        append!(traj, Trajectory(dm, X0[i], T1[i], T2[i]))
    end

    return traj
end

"""
    copy(traj::Trajectory)

Return a shallow copy of traj
"""
function Base.copy(traj::Trajectory)
    # Trajectory(dm::DynamicalModel, X0, T)

    x0s = Vector{Vector}()
    ts = Vector{Vector}()
    for arc in traj
        push!(x0s, arc[begin])
        push!(ts, [arc.t[begin], arc.t[end]])
    end

    return Trajectory(dm(traj), x0s, ts)
end


"""
    dimension(traj::Trajectory{D}) where D

Return dimension of the dynamical model of traj
"""
dimension(traj::Trajectory{D}) where {D} = D

"""
    x0(traj::Trajectory)

Return the initial condition of the trajectory
"""
x0(traj::Trajectory) = traj.X_0

"""
    dm(traj::Trajectory)

Return the dynamical model of the trajectory
"""
dm(traj::Trajectory) = traj.dm

"""
    tspan(traj::Trajectory)

Return the time span of the trajectory
"""
tspan(traj::Trajectory) = traj.tspan
tspan(sol::OrdinaryDiffEq.ODESolution) = sort([sol.t[begin], sol.t[end]])

"""
    tof(traj::Trajectory)

Return the time of flight of the trajectory
"""
tof(traj::Trajectory) = traj.tspan[2] - traj.tspan[1]
tof(sol::OrdinaryDiffEq.ODESolution) = tspan(sol)[end] - tspan(sol)[begin]

"""
    solvec(traj::Trajectory)

Return the vector of ODESolutions that defines the trajectory
"""
solvec(traj::Trajectory) = traj.X

"""
    length(traj::Trajectory)

Return the number of ODESolutions comprising the trajectory
"""
Base.length(traj::Trajectory) = length(solvec(traj))

"""
    getindex(traj::Trajectory, i::Int)

Return i'th element in X
"""
Base.getindex(traj::Trajectory, i::Int) = traj.X[i]

function Base.getindex(traj::Trajectory, r::UnitRange{Int})
    outvec = Vector{eltype(solvec(traj))}()
    for i in r
        push!(outvec, traj[i])
    end
    return outvec
end


"""
    iterate(::Trajectory)

Method for iterating through FreeVariables
"""
Base.iterate(traj::Trajectory, state=1) = state>length(traj) ? nothing : (traj[state], state+1)

"""
    Base.append!(traj1::Trajectory, traj2::Trajectory)

Extend Base.append! for trajectories.
"""
function Base.append!(traj1::Trajectory, trajn...)
    # Loop through input trajectories
    for i = 1:length(trajn)
        # Check that trajn consists of Trajectory objects with the same dynamical model
        if typeof(trajn[i]) <: Trajectory && dm(traj1) == dm(trajn[i])
            # Recompute trajectory if initial time of trajn[i] is not 
            # the same as the final time of traj1
            if tspan(traj1)[end] != tspan(trajn[i])[1]
                # Preserve ODESolution segments in trajn[i]

                # tstart = tspan(traj1)[end]
                tstart = maximum(tspan(traj1))
                for j = 1:length(trajn[i])
                    # traj2 = Trajectory(dm(trajn[i]), trajn[i].X[j][begin], [0.0, tof(trajn[i].X[j])].+tstart)
                    traj2 = Trajectory(dm(trajn[i]), get_x0(trajn[i][j]), [0.0, tof(trajn[i].X[j])].+tstart)
                    append!(traj1.X, traj2.X)
                    tstart = tstart + tof(trajn[i].X[j])
                end

                # println("reprop")
                # traj2 = Trajectory(dm(trajn[i]), x0(trajn[i]), [0.0, tof(trajn[i])].+tspan(traj1)[end])
                # append!(traj1.X, traj2.X)
            else
                append!(traj1.X, trajn[i].X)
            end


            traj1.tspan[2]+=tof(trajn[i])
        else
            throw(TypeError(:append!, "", Trajectory, typeof(trajn[i])))
        end
    end

end

"""
    __within(x::T, left::T, right::T)

Check if x is within the range between left and right
"""
function __within(x::Real, left::Real, right::Real)
    # Check that left < right
    if left > right
        # Swap left and right
        temp = left
        left = right
        right = temp # old left
    end
    if x <= right && x >= left
        return true
    else
        return false
    end
end

"""
    (traj::Trajectory)(T::Real)

Return the state along the trajectory at time T or times "times"
"""
function (traj::Trajectory)(T::Real)
    # Check that T is within tspan
    if !__within(T, tspan(traj)[1], tspan(traj)[2])
        if isapprox(tspan(traj)[1],T) 
            T = tspan(traj)[1]

        elseif isapprox(tspan(traj)[2],T)
            T = tspan(traj)[2]
        else
            throw(ErrorException("T = $(T) is not within tspan, $(tspan(traj))"))
        end
    end
    
    # Return the state at time T
    for i = 1:length(traj)
        if __within(T, solvec(traj)[i].t[begin], solvec(traj)[i].t[end]+DEFAULT_CONVERGENCE_TOL)
            return solvec(traj)[i](T)
        end
    end

    # Allow for slight numerical errors in individual solvec tspans
    # and output final time of trajectory even if it is technically
    # out of range by the above loop
    if abs(T - tspan(traj)[2]) < DEFAULT_CONVERGENCE_TOL
        return solvec(traj)[end].u[end]
    end
    if abs(T - tspan(traj)[1]) < DEFAULT_CONVERGENCE_TOL
        return solvec(traj)[begin].u[begin]
    end

    # TODO fix case where a requested time lies between two segments -- 
    # TODO check that fix on line 256 works

    throw(ErrorException("T is not within any element of traj.X"))
end

function (traj::Trajectory)(times) 
    # outvec = Vector{Vector{Real}}()
    outvec = Vector{Vector{Float64}}()
    foreach(x->push!(outvec, traj(x)), times)
    return outvec
end

"""
    isperiodic(traj::Trajectory, tol=DEFAULT_CONVERGENCE_TOL)

Returns true if the final state and beginning state in traj are the same, within a tolerance
"""
function isperiodic(traj::Trajectory, tol=DEFAULT_CONVERGENCE_TOL)
    abserr = map(x->abs(x),traj(tspan(traj)[1]) - traj(tspan(traj)[2]))

    if all(abserr.<tol)
        return true
    else
        return false
    end
end

"""
    iscontinuous(traj::Trajectory, tol=DEFAULT_CONVERGENCE_TOL)

Check that all segments in traj flow into the next one continuously, with a tolerance
"""
function iscontinuous(traj::Trajectory, tol=DEFAULT_CONVERGENCE_TOL; inds=:)
    if length(traj) == 1
        return true
    end
    for i = 1:length(traj)-1
        # abserr = map(x->abs(x),traj[i+1][begin]-traj[i][end])
        abserr = map(x->abs(x),traj[i+1](minimum(traj[i+1].t))[inds]-traj[i](maximum(traj[i].t))[inds])
        if all(abserr.<tol)
            # return true
            continue
        else
            println("Discontinuous between segments $(i) and $(i+1)")
            return false
        end

    end

    return true
end


"""
    get_x0(sol::OrdinaryDiffEq.ODESolution)

Return the state corresponding to the smallest time
in sol.t
"""
function get_x0(sol::OrdinaryDiffEq.ODESolution)
    if sol.t[begin] < sol.t[end]
        return sol.u[begin]
    else
        return sol.u[end]
    end
end

function get_x0(solvec::Vector{OrdinaryDiffEq.ODESolution})
    outvec = []
    for i = 1:length(solvec)
        push!(outvec, get_x0(solvec[i]))
    end

    return outvec
end
get_x0(traj::Trajectory) = get_x0(solvec(traj))

"""
    get_xf(sol::OrdinaryDiffEq.ODESolution)

Return the state corresponding to the largest time
in sol.t
"""
function get_xf(sol::OrdinaryDiffEq.ODESolution)
    if sol.t[begin] < sol.t[end]
        return sol.u[end]
    else
        return sol.u[begin]
    end
end
function get_xf(solvec::Vector{OrdinaryDiffEq.ODESolution})
    outvec = []
    for i = 1:length(solvec)
        push!(outvec, get_xf(solvec[i]))
    end

    return outvec
end
get_xf(traj::Trajectory) = get_xf(solvec(traj))


"""
    get_t0(sol::OrdinaryDiffEq.ODESolution)

Return the smallest time in sol
"""
function get_t0(sol::OrdinaryDiffEq.ODESolution)
    if sol.t[begin] < sol.t[end]
        return sol.t[begin]
    else
        return sol.t[end]
    end
end
function get_t0(solvec::Vector{OrdinaryDiffEq.ODESolution})
    outvec = []
    for i = 1:length(solvec)
        push!(outvec, get_t0(solvec[i]))
    end

    return outvec
end
get_t0(traj::Trajectory) = get_t0(solvec(traj))

"""
    get_tf(sol::OrdinaryDiffEq.ODESolution)

Return the largest time in sol
"""
function get_tf(sol::OrdinaryDiffEq.ODESolution)
    if sol.t[begin] < sol.t[end]
        return sol.t[end]
    else
        return sol.t[begin]
    end
end
function get_tf(solvec::Vector{OrdinaryDiffEq.ODESolution})
    outvec = []
    for i = 1:length(solvec)
        push!(outvec, get_tf(solvec[i]))
    end

    return outvec
end
get_tf(traj::Trajectory) = get_tf(solvec(traj))


"""
    stm(traj::Trajectory{D}) where {D}

Return the state transition matrix of the trajectory, from the initial
time to a final time (which is the TOF of the trajectory by default)
"""
function stm(traj::Trajectory{D}, T::Real=tof(traj)) where {D}
    mats = Vector{AbstractMatrix}()

    for arc in traj
        if __within(T, arc.t[begin], arc.t[end])
            sol = tangent_solve(dm(traj), arc[begin], (arc.t[begin], T))
            # println(size(sol.u[end]))
            push!(mats, sol.u[end][:,2:end])
            break # TODO need a more elegant solution
        else
            sol = tangent_solve(dm(traj), arc[begin], (arc.t[begin], arc.t[end]))
            # println(size(sol.u[end]))
            push!(mats, sol.u[end][:,2:end])
        end
    end

    Φ = Matrix{Float64}(mats[end])
    for i = length(mats)-1:-1:1
        Φ = Φ*mats[i]
    end
    
    return Φ

end

"""
    Base.show

Overload the show operator to pretty print the Trajectory to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", traj::Trajectory{D}) where {D}
    print(io, "Trajectory\n")
    print(io, "- Dimension: $(D)\n")
    print(io, "- Length: $(length(traj))\n")
    print(io, "- Time Span: $(tspan(traj))\n")
    print(io, "- X0: $(x0(traj))\n")
end

"""
    Base.show

Overload the show operator to pretty print the Vector{ODESolution} to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", solvec::Vector{OrdinaryDiffEq.ODESolution})
    for i = 1:length(solvec)
        println("\n\nODESolution Number: $(i)\n≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡")
        if isassigned(solvec,i)
            show(io, "text/plain", solvec[i])
        else
            println("#undef")
        end
    end
end

"""
    to_dict(traj::Trajectory)

Convert Trajectory `traj` into a Dict
"""
function to_dict(traj::Trajectory; continuity_tol=DEFAULT_CONVERGENCE_TOL)
    model = dm(traj) # Dynamical model of the periodic orbit

    outdict = Dict() # Dictionary that will be used to output results

    outdict["datatype"] = "Trajectory"
    outdict["iscontinuous"] = iscontinuous(traj, continuity_tol)


    # Add model information to the dictionary
    if typeof(model) <: Cr3bpModel

        ## Write CR3BP model info
        outdict["modeltype"] = "Cr3bpModel"
        if isnothing(model.primaries)
            outdict["P1"] = ""
            outdict["P2"] = ""
            outdict["mu"] = mass_ratio(model)

        else
            outdict["P1"] = primary_bodies(model)[1].name
            outdict["P2"] = primary_bodies(model)[2].name
            outdict["mu"] = mass_ratio(model)
        end
        
    elseif typeof(model) <: HFEModel
        epoch = get_epoch(model)
        outdict["modeltype"] = "HFEModel"
        outdict["central_body"] = get_central_body(model).name
        outdict["additional_bodies"] = map(x->x.name, get_additional_bodies(model))

        outdict["epoch_type"] = string(typeof(epoch))
        outdict["epoch_year"] = epoch.year
        outdict["epoch_month"] = epoch.month
        outdict["epoch_day"] = epoch.day
        outdict["epoch_hour"] = epoch.hour
        outdict["epoch_minute"] = epoch.minute
        outdict["epoch_second"] = epoch.second

        outdict["dimensional_mass"] = dimensional_mass(model)
        outdict["dimensional_length"] = dimensional_length(model)
        outdict["dimensional_time"] = dimensional_time(model)
    end

    ## Get states and times
    outdict["X0"] = Vector{Vector{Float64}}(get_x0(traj))
    outdict["T0"] = Vector{Float64}(get_t0(traj))
    outdict["Tf"] = Vector{Float64}(get_tf(traj))

    return outdict
end

function to_dict(trajvec::Vector{Trajectory})
    outdict = Dict("data"=>[])
    for traj in trajvec
        push!(outdict["data"], to_dict(traj))
    end
    
    return outdict
end

"""
    to_mat(traj::Trajectory, filename::String)

Save Trajectory `traj` to a .mat file 
"""
function to_mat(traj::Trajectory, filename::String)
    outdict = to_dict(traj)
    endswith(filename, ".mat") ? nothing : filename = filename*".mat"
    matwrite(filename,Dict("data"=>[outdict]))
    println("Saved to $(filename)")
end

function to_mat(trajvec::Vector{Trajectory}, filename::String)
    outdict = to_dict(trajvec)
    endswith(filename, ".mat") ? nothing : filename = filename*".mat"
    matwrite(filename,outdict)
    println("Saved to $(filename)")
end
