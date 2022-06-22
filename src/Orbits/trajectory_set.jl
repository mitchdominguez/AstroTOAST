# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                 TRAJECTORY SET
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
    struct TrajectorySet{D,N}

Object that represents a collection of trajectories that all start at the same epoch
and have the same time of flight
"""

struct TrajectorySet{D,N}
    X_0::Vector{Float64}
    tspan::Vector{Float64}
    dm::DynamicalModel 
    trajvec::Vector{Trajectory{D}}

    # Constructor for generating a trajectory from a single patch point
    function TrajectorySet(dm::DynamicalModel, trajvec::Vector{Trajectory{D}}) where {D}
        # Check that dimensions match
        dim = dimension(dm)
        if D != dim
            throw(DimensionMismatch("trajectory and dynamical model have different dimension"))
        end

        # Check that time of flight matches
        # and build X_0
        proptimes = tspan(trajvec[1])
        X_0 = Vector{Float64}()
        for traj in trajvec
            if tspan(traj) != proptimes
                throw(ErrorException("All trajectories must have the same time of flight"))
            end
            append!(X_0, x0(traj))
        end

        N = length(trajvec)
        new{D,N}(X_0, proptimes, dm, trajvec)
    end
end

########################## OUTER CONSTRUCTORS ##########################
"""
    TrajectorySet(dm::DynamicalModel, X0::FreeVariable, T::FreeVariable)

Take in a dynamical model, a vector `X0` which holds `N` initial conditions,
and generates a TrajectorySet with `T` time of flight
"""
function TrajectorySet(dm::DynamicalModel, X0::Vector{Float64}, T::Float64; backprop=false)
    D = dimension(dm)

    # Check that there are an integer number of initial condition sets within X0
    if length(X0) % D != 0
        throw(DimensionMismatch("X0 must be of full_length $(D)N"))
    end

    N = Int(length(X0)/D)

    trajvec = Vector{Trajectory{D}}(undef,N)

    # Create vector of trajectories
    for i = 1:N
        q0 = X0[D*i-(D-1):D*i]

        if backprop && typeof(dm)<:Cr3bpModel
            # Backwards propagation
            trajvec[i] = Trajectory(dm, q0, -T)

        elseif !backprop && typeof(dm)<:Cr3bpModel
            # Forwards propagation
            trajvec[i] = Trajectory(dm, q0, T)
        end

    end

    return TrajectorySet(dm, trajvec)
end

"""
    dimension(ts::TrajectorySet{D,N}) where {D,N}

Return dimension of the dynamical model of traj
"""
dimension(ts::TrajectorySet{D,N}) where {D,N} = D

"""
    numels(ts::TrajectorySet{D,N}) where {D,N}

Return number of trajectories in `ts`
"""
numels(ts::TrajectorySet{D,N}) where {D,N} = N
Base.length(ts::TrajectorySet{D,N}) where {D,N} = N

"""
    x0(ts::TrajectorySet)

Return the initial condition of the trajectory set
"""
x0(ts::TrajectorySet) = ts.X_0

"""
    dm(ts::TrajectorySet)

Return the dynamical model of the trajectory set
"""
dm(ts::TrajectorySet) = ts.dm

"""
    tspan(ts::TrajectorySet)

Return the time span of the trajectory set
"""
tspan(ts::TrajectorySet) = ts.tspan

"""
    tof(ts::TrajectorySet)

Return the time of flight of the trajectory set
"""
tof(ts::TrajectorySet) = ts.tspan[2] - ts.tspan[1]

"""
    getindex(ts::TrajectorySet, i::Int)

Return i'th element in X
"""
Base.getindex(ts::TrajectorySet, i::Int) = ts.trajvec[i]

function Base.getindex(ts::TrajectorySet, r::UnitRange{Int})
    outvec = Vector{Trajectory}()
    for i in r
        push!(outvec, ts[i])
    end
    return outvec
end

"""
    iterate(ts::TrajectorySet)

Method for iterating through Trajectories in a TrajectorySet
"""
Base.iterate(ts::TrajectorySet, state=1) = state>numels(ts) ? nothing : (ts[state], state+1)

"""
    (ts::TrajectorySet)(T::Real)

Return the states along the trajectory at time T or times "times"
"""
function (ts::TrajectorySet{D,N})(T::Real) where {D,N}
    # Check that T is within tspan
        if isapprox(tspan(ts)[1],T) 
            T = tspan(ts)[1]

        elseif isapprox(tspan(ts)[2],T)
            T = tspan(ts)[2]
        else
            throw(ErrorException("T = $(T) is not within tspan, $(tspan(ts))"))
        end
    
    # Return the state at time T
    xT = similar(x0(ts))
    for i = 1:N
        xT[D*i-(D-1):D*i] = ts[i](T)
    end

    return xT

    throw(ErrorException("T was not within any element of traj.X"))
end

function (ts::TrajectorySet)(times) 
    outvec = Vector{Vector{Real}}()
    foreach(x->push!(outvec, ts(x)), times)
    return outvec
end

"""
    Base.show

Overload the show operator to pretty print the PeriodicOrbit to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", ts::TrajectorySet{D,N}) where {D,N}
    print(io, "Trajectory Set: \n")
    print(io, "- Dimension: $(D)\n")
    print(io, "- Number of Trajectories: $(N)\n")
    print(io, "- Time Span: $(ts.tspan) \n")
end
