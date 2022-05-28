# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                    PERIODIC ORBIT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
using LinearAlgebra

"""
    struct PeriodicOrbit{D}

Object that represents periodic trajectories within a dynamical model
"""
struct PeriodicOrbit{D}
    traj::Trajectory{D}
    M::Matrix{Float64}
    λ::Vector{ComplexF64}
    V::Vector{Vector{ComplexF64}}
    name::String
    family::String

    function PeriodicOrbit(traj::Trajectory{D}, name = "", family = "", tol=DEFAULT_ABS_TOL) where {D}
        # Ensure that traj is periodic
        if !isperiodic(traj)
            throw(ErrorException("traj is not periodic!"))
        end

        # Ensure that all segments of traj are continuous
        if !iscontinuous(traj)
            throw(ErrorException("traj is not continuous!"))
        end

        # Calculate monodromy matrix
        M = stm(traj)
        
        # Calculate eigenvalues and eigenvectors
        λ, vec = eigen(M)

        # Check that eigenvalues and eigenvectors are paired properly
        err = Vector{Float64}(undef, D)
        for i = 1:D
            err[i] = norm((M-I(D)*λ[i])*vec[:,i])
        end
        if !all(err.<tol)
            throw(ErrorException("Eigenvalues not paired properly"))
            # TODO make a function to pair eigenvalues instead of just
            # erroring out
        end

        V = Vector{Vector{ComplexF64}}(undef,D)
        for i = 1:D
            V[i] = vec[:,i]
        end

        # Create new PeriodicOrbit
        new{D}(copy(traj), M, λ, V, name, family)

    end
end

############################
# General Utilities
############################
"""
    traj(po::PeriodicOrbit)

Return the Trajectory of the PO
"""
traj(po::PeriodicOrbit) = po.traj

"""
    dimension(po::PeriodicOrbit{D}) where D

Return dimension of the dynamical model of traj
"""
dimension(po::PeriodicOrbit{D}) where {D} = D

"""
    x0(po::PeriodicOrbit)

Return the initial condition of the trajectory
"""
x0(po::PeriodicOrbit) = x0(traj(po))

"""
    dm(po::PeriodicOrbit)

Return the dynamical model of the trajectory
"""
dm(po::PeriodicOrbit) = dm(traj(po))

"""
    period(po::PeriodicOrbit)

Return the time of flight of the trajectory
"""
period(po::PeriodicOrbit) = tof(traj(po))

"""
    jacobi_constant(po::PeriodicOrbit)

Return the jacobi constant of the periodic orbit
"""
jacobi_constant(po::PeriodicOrbit) = typeof(dm(po)) <: Cr3bpModel ? jacobi_constant(dm(po), x0(po)) : throw(MethodError(jacobi_constant, po))

"""
    (po::PeriodicOrbit)(T)

Return the state on the periodic orbit at time T
"""
(po::PeriodicOrbit)(T) = traj(po)(T)

"""
    monodromy(po::PeriodicOrbit)

Return the monodromy matrix of the periodic orbit
"""
monodromy(po::PeriodicOrbit) = po.M

"""
    name(po::PeriodicOrbit)

Print the name of the Periodic Orbit
"""
name(po::PeriodicOrbit) = po.name

"""
    family(po::PeriodicOrbit)

Print the family of the Periodic Orbit
"""
family(po::PeriodicOrbit) = po.family

######################################################
# Stability, eigendecomoposition, manifolds
######################################################

"""
    λ(po::PeriodicOrbit)
"""
#TODO make sure 9:2 NRHO has the correct synodic period
#TODO periapsis, apoapsis
#TODO eigenvalue/vector sorting
#TODO stable/unstable/center eigenvalue/vector output

"""
    Base.show

Overload the show operator to pretty print the PeriodicOrbit to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", po::PeriodicOrbit{D}) where {D}
    print(io, "Periodic Orbit: $(name(po)) ($(family(po)))\n")
    print(io, "- Dimension: $(D)\n")
    print(io, "- Period: $(period(po)) ndim = $(period(po)*dimensional_time(dm(po))*sec2day) days\n")
    print(io, "- X0: $(x0(po))\n")
end
