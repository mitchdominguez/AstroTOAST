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
    JC::Float64

    function PeriodicOrbit(traj::Trajectory{D}, tol=DEFAULT_ABS_TOL) where {D}
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
        end

        V = Vector{Vector{ComplexF64}}(undef,D)
        for i = 1:D
            V[i] = vec[:,i]
        end

        # Calculate Jacobi constant
        JC = jacobi_constant(dm(traj), x0(traj))

        # Create new PeriodicOrbit
        new{D}(traj, M, λ, V, JC)

    end
end

# """
    # Base.show

# Overload the show operator to pretty print the Trajectory to the console.
# """
# function Base.show(io::IO, ::MIME"text/plain", traj::Trajectory{D}) where {D}
    # print(io, "Trajectory\n")
    # print(io, "- Dimension: $(D)\n")
    # print(io, "- Length: $(length(traj))\n")
    # print(io, "- Time Span: $(traj.tspan)\n")
    # print(io, "- X0: $(traj.X_0)\n")
# end
