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
        sol = solve(dm, X0, tspan)

        new{dimension(dm)}(X0, tspan, dm, [sol])
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
            tspan(trajn[i])
            if tspan(traj1)[end] != tspan(trajn[i])[1]
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
    within(x::T, left::T, right::T)

Check if x is within the range between left and right
"""
function __within(x::Real, left::Real, right::Real)
    # Check that left < right
    if left > right
        throw(ErrorException("left must be <= right"))
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
        throw(ErrorException("T is not within tspan"))
    end
    
    # Return the state at time T
    for i = 1:length(traj)
        if __within(T, solvec(traj)[i].t[begin], solvec(traj)[i].t[end])
            return solvec(traj)[i](T)
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
    isperiodic(traj::Trajectory, tol=1e-12)

Returns true if the final state and beginning state in traj are the same, within a tolerance
"""
function isperiodic(traj::Trajectory, tol=DEFAULT_ABS_TOL)
    abserr = map(x->abs(x),traj(tspan(traj)[1]) - traj(tspan(traj)[2]))

    if all(abserr.<tol)
        return true
    else
        return false
    end
end


"""
    Base.show

Overload the show operator to pretty print the Trajectory to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", traj::Trajectory{D}) where {D}
    print(io, "Trajectory\n")
    print(io, "- Dimension: $(D)\n")
    print(io, "- Length: $(length(traj))\n")
    print(io, "- Time Span: $(traj.tspan)\n")
    print(io, "- X0: $(traj.X_0)\n")
end

"""
    Base.show

Overload the show operator to pretty print the Vector{ODESolution} to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", solvec::Vector{OrdinaryDiffEq.ODESolution})
    for i = 1:length(solvec)
        println("ODESolution Number: $(i)\n≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡")
        if isassigned(solvec,i)
            show(io, "text/plain", solvec[i])
        else
            println("#undef")
        end
    end
end
