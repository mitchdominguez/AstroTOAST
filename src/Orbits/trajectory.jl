# TODO export initial conditions for each segment with TOFs



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
    if length(X0) != length(T)
        throw(DimensionMismatch("Different numbers of patch points and times given"))
    end

    len = length(X0)

    traj = Trajectory(dm, X0[1], T[1])

    for i = 2:len
        append!(traj, Trajectory(dm, X0[i], T[i]))
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

"""
    tof(traj::Trajectory)

Return the time of flight of the trajectory
"""
tof(traj::Trajectory) = traj.tspan[2] - traj.tspan[1]

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
    # TODO make append! work for trajectories with multiple segments
    # Loop through input trajectories
    for i = 1:length(trajn)
        # Check that trajn consists of Trajectory objects with the same dynamical model
        if typeof(trajn[i]) <: Trajectory && dm(traj1) == dm(trajn[i])
            # Recompute trajectory if initial time of trajn[i] is not 
            # the same as the final time of traj1
            tspan(trajn[i])
            if tspan(traj1)[end] != tspan(trajn[i])[1]
                # println("reprop")
                traj2 = Trajectory(dm(trajn[i]), x0(trajn[i]), [0.0, tof(trajn[i])].+tspan(traj1)[end])
                append!(traj1.X, traj2.X)
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
        if __within(T, solvec(traj)[i].t[begin], solvec(traj)[i].t[end])
            return Vector(solvec(traj)[i](T))
        end
    end

    throw(ErrorException("T was not within any element of traj.X"))
end

function (traj::Trajectory)(times) 
    outvec = Vector{Vector{Real}}()
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
function iscontinuous(traj::Trajectory, tol=DEFAULT_CONVERGENCE_TOL)
    if length(traj) == 1
        return true
    end
    for i = 1:length(traj)-1
        # abserr = map(x->abs(x),traj[i+1][begin]-traj[i][end])
        abserr = map(x->abs(x),traj[i+1](minimum(traj[i+1].t))-traj[i](maximum(traj[i].t)))
        if all(abserr.<tol)
            return true
        else
            return false
        end

    end

end


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
