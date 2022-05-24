# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                     TRAJECTORY
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
# using OrdinaryDiffEq

"""
    struct Trajectory{D}

Object that will represent trajectories within a given dynamical model
"""

struct Trajectory{D}
    X_0::Vector{Real}
    tspan::NTuple{2,Real}
    dm::DynamicalModel 
    X::Vector{OrdinaryDiffEq.ODESolution} # State history

    function Trajectory(dm::DynamicalModel, X0::AbstractVector{T}, proptime::Union{T, Vector{T}, NTuple{2,Real}}) where {T<:Real}
        # Check that dimensions all match
        if dimension(dm) != length(X0)
            throw(DimensionMismatch("Initial conditions and dynamical model have different dimensions"))
        end
        
        # Create tspan 
        if length(proptime) == 1
            tspan = promote(T(0.0), proptime[1])
        elseif length(proptime) == 2
            tspan = promote(T(proptime[1]), proptime[2])
        else
            throw(DimensionMismatch("proptime must be of length 1 or 2"))
        end

        # Propagate initial conditions for the desired time
        sol = solve(dm, X0, tspan)

        new{dimension(dm)}(X0, tspan, dm, [sol])
    end
end

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
    within(x::T, left::T, right::T)

Check if x is within the range between left and right
"""
function within(x::Real, left::Real, right::Real)
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

Return the state along the trajectory at time T
"""
function (traj::Trajectory)(T::Real)
    # Check that T is within tspan
    if !within(T, tspan(traj)[1], tspan(traj)[2])
        throw(ErrorException("T is not within tspan"))
    end
    
    # Return the state at time T
    for i = 1:length(traj)
        if within(T, solvec(traj)[i].t[begin], solvec(traj)[i].t[end])
            return solvec(traj)[i](T)
        end
    end

    throw(ErrorException("T was not within any element of traj.X"))
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
        show(io, "text/plain", solvec[i])
    end
end
