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
    traj::Trajectory{D} # Trajectory
    M::Matrix{Float64} # Monodromy matrix
    λ::Vector{ComplexF64} # Eigenvalues
    V::Vector{Vector{ComplexF64}} # Eigenvectors
    name::String # Name of orbit
    family::String # Family that orbit belongs to

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

        V = Vector{Vector{ComplexF64}}(undef,D)
        for i = 1:D
            V[i] = vec[:,i]
        end

        # Make sure the eigenvalues and eigenvectors are paired and sorted properly
        sorteigs!(V, λ, M)

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
    eigvals(po::PeriodicOrbit)

Return eigenvalues of the PeriodicOrbit
"""
LinearAlgebra.eigvals(po::PeriodicOrbit) = copy(po.λ)

"""
    eigvecs(po::PeriodicOrbit)

Return eigenvectors of the PeriodicOrbit
"""
LinearAlgebra.eigvecs(po::PeriodicOrbit) = copy(po.V)

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
    paireigs(v::Vector{Vector{ComplexF64}}, M::Matrix{Float64}, λ::Vector{ComplexF64})

Pair eigenvalues with their corresponding eigenvectors.

This function will change the value of v in place
"""
function paireigs!(v::Vector{Vector{ComplexF64}}, λ::Vector{ComplexF64}, M::Matrix{Float64})
    vcopy = copy(v)
    D = length(λ)

    err = Vector{Float64}(undef, D)

    for i = 1:D
        # Calculate error for each eigenvector
        for j = 1:length(vcopy)
            err[j] = norm((M-I(D)*λ[i])*vcopy[j])
        end
        
        # Find element of vcopy that results in the minimum error
        minval, ind = findmin(err)

        # Write eigenvector corresponding to minval to index i in v
        v[i] = vcopy[ind]

        # Remove the element from vcopy and err
        deleteat!(err, ind)
        deleteat!(vcopy, ind)
    end

    return v
end

"""
    sorteigs(v::Vector{Vector{ComplexF64}}, λ::Vector{ComplexF64}, M::Matrix{Float64})

Sort pairs of eigenvalues from the pairs with the highest magnitude to pairs with the magnitude closest to 1

This function will change the value of λ and v in place
"""
function sorteigs!(v::Vector{Vector{ComplexF64}}, λ::Vector{ComplexF64}, M::Matrix{Float64})
    # Error if odd number of eigenvalues
    if length(λ)%2 != 0
        throw(DimensionMismatch("Odd number of eigenvalues"))
    end

    # Pairs of eigenvalues
    numpairs = Int(length(λ)/2)

    # Copy λ
    lamcopy = copy(λ)


    for i = 1:2:numpairs*2
        # Sort by eigenvalue magnitude
        mags = map(x->abs(x), lamcopy)

        # Find maximum 
        maxval, maxind = findmax(mags)

        # Find reciprocal of lamcopy[ind]
        recip, recipind = findmin(map(x->abs(x), lamcopy.^(-1) .- lamcopy[maxind]))

        # Write lamcopy[maxind] and lamcopy[recipind] into elements i, i+1 in λ
        λ[i] = lamcopy[maxind]
        λ[i+1] = lamcopy[recipind]

        # Remove lamcopy[maxind] and lamcopy[recipind] from lamcopy
        deleteat!(lamcopy, Tuple([x for x in sort([maxind,recipind])]))

    end

    # Re-pair eigenvalues and eigenvectors
    paireigs!(v, λ, M)
    

    return (λ, v)
end


#TODO periapsis, apoapsis
#TODO stable/unstable/center eigenvalue/vector output
#TODO number of stable/unstable/center eigs

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
