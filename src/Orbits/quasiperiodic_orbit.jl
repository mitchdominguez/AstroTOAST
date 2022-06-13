# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                               QUASIPERIODIC ORBIT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
using LinearAlgebra

"""
    struct QuasiPeriodicOrbit{D}

Object that represents periodic trajectories within a dynamical model
"""
struct QuasiPeriodicOrbit{D,N}
    ic::TrajectorySet{D,N} # Invariant curve
    T::Float64 # Stroboscopic time
    ρ::Float64 # Twist angle
    fixedpt::Vector{Float64} # Fixed point used to generate invariant curve # TODO make this a PO
    DG::Matrix{Float64} # DG matrix
    λ::Vector{ComplexF64} # Eigenvalues of DG
    V::Vector{Vector{ComplexF64}} # Eigenvectors of DG
    name::String # Name of QPO
    family::String # Family that QPO belongs to

    function QuasiPeriodicOrbit(ts::TrajectorySet{D,N}, rho::Float64, xstar::Vector{Float64}, 
            name = "", family = "", tol=DEFAULT_CONVERGENCE_TOL) where {D,N}
        # Check that fixed point has the right dimension
        if length(xstar) != D
            throw(DimensionMismatch("dimension of xstar, and DynamicalModel are incompatible"))
        end

        # Check that N is odd
        if iseven(N)
            throw(InvalidStateException("There must be an odd number of nodes in the invariant curve", :N))
        end

        # Check that rho and T are both 1 dimensional
        if length(rho) != 1
            throw(DimensionMismatch("rho must have full length 1"))
        end
        
        ############### Check that each point on the invariant curve has the same jacobi constant ###############
        model = dm(ts)
        jc = jacobi_constant(model, x0(ts[1]))
        for traj in ts
            if abs(jacobi_constant(model, x0(traj)) - jc) > tol
                throw(ErrorException("Incompatible Jacobi constants on invariant curve"))
            end
        end
        

        ############### Check that the invariance constraint is met ###############
        # Generate FreeVariable for invariant curve
        x0vec = x0(ts)
        u0vec = similar(x0(ts))
        for i = 1:N
            u0vec[D*i-(D-1):D*i] = x0vec[D*i-(D-1):D*i] - xstar # Full state including fixed point
        end
        U0 = FreeVariable("U0", u0vec)

        # Generate FreeVariable for stroboscopic time
        T = FreeVariable("T", tof(ts))

        # Generate FreeVariable for twist angle
        rhovar = FreeVariable("rho", rho)

        # Invariance Constraint
        ic = InvarianceConstraint2D(U0, xstar, T, rhovar, dm(ts))

        if norm(evalconstraint(ic)) > tol
            println(norm(evalconstraint(ic)))
            throw(ErrorException("The given states do not correspond to a 2D torus, given the tolerance provided"))
        end

        # Generate DG matrix and its eigendecomposition
        DG = __dIC_du0{D*N}()(ic) + I(D*N)

        lam, vee = eigen(DG)

        V = Vector{Vector{ComplexF64}}(undef,D*N)
        for i = 1:D*N
            V[i] = vee[:,i]
        end

        return new{D,N}(ts, tof(ts), rho, xstar, DG, lam, V, name, family)
    end
end

########################## OUTER CONSTRUCTORS ##########################
function QuasiPeriodicOrbit(dm::DynamicalModel, X0::Vector{Float64}, T::Float64,  
        rho::Float64, xstar::Vector{Float64}; 
        name = "", family = "", tol=DEFAULT_CONVERGENCE_TOL) 
    ts = TrajectorySet(dm, X0, T)
    return QuasiPeriodicOrbit(ts, rho, xstar, name, family, tol)
end

########################## Basic Utilities ##########################

"""
    x0(qpo::QuasiPeriodicOrbit)

Return the initial states on the invariant curve as a Vector
"""
x0(qpo::QuasiPeriodicOrbit) = x0(qpo.ic)

"""
    xstar(qpo::QuasiPeriodicOrbit)

Return the fixed point at the initial states of the invariant curve of the QPO
"""
xstar(qpo::QuasiPeriodicOrbit) = qpo.fixedpt

"""
    u0(qpo::QuasiPeriodicOrbit)

Return the initial states on the invariant curve relative to the fixed point
"""
function u0(qpo::QuasiPeriodicOrbit{D,N}) where {D,N} 
    x0vec = x0(qpo)
    u0vec = similar(x0vec)
    for i = 1:N
        u0vec[D*i-(D-1):D*i] = x0vec[D*i-(D-1):D*i] - xstar(qpo) # Full state including fixed point
    end
    return u0vec
end

"""
    u2x(u::AbstractVector, xfixed::AbstractVector)

Add fixed point `xfixed` to a relative state on the qpo `u`
"""
function u2x(u::AbstractVector, xfixed::AbstractVector)
    D = length(xfixed)

    if length(u)%D != 0
        throw(ErrorException("Lengths of u and xfixed are incompatible"))
    end

    N = Int(length(u)/D)

    xvec = similar(u)
    for i = 1:N
        xvec[D*i-(D-1):D*i] = u[D*i-(D-1):D*i] + xfixed # Full state including fixed point
    end
    return xvec

end

"""
    x2u(x::AbstractVector, xfixed::AbstractVector)

Add fixed point `xfixed` to a relative state on the qpo `u`
"""
function x2u(x::AbstractVector, xfixed::AbstractVector)
    D = length(xfixed)

    if length(x)%D != 0
        throw(ErrorException("Lengths of u and xfixed are incompatible"))
    end

    N = Int(length(x)/D)

    uvec = similar(x)
    for i = 1:N
        uvec[D*i-(D-1):D*i] = x[D*i-(D-1):D*i] - xfixed # Full state including fixed point
    end
    return uvec

end

"""
    strobetime(qpo::QuasiPeriodicOrbit)

Return the stroboscopic period of the torus to be targeted
"""
strobetime(qpo::QuasiPeriodicOrbit) = qpo.T

"""
    rotationangle(qpo::QuasiPeriodicOrbit)

Return the angle by which the invariant curve will rotate after one period T
"""
rotationangle(qpo::QuasiPeriodicOrbit) = qpo.ρ

"""
    invariantcurve(qpo::QuasiPeriodicOrbit)

Return the TrajectorySet representing the invariant curve of the QPO
"""
invariantcurve(qpo::QuasiPeriodicOrbit) = qpo.ic

"""
    dm(qpo::QuasiPeriodicOrbit)

Return the dynamical model of the trajectory
"""
dm(qpo::QuasiPeriodicOrbit) = dm(invariantcurve(qpo))

"""
    dimension(qpo::QuasiPeriodicOrbit)

Return the dimension of the dynamical model of the QPO
"""
dimension(qpo::QuasiPeriodicOrbit{D,N}) where {D,N} = D

"""
    numels(qpo::QuasiPeriodicOrbit)

Return the dimension of the dynamical model of the QPO
"""
numels(qpo::QuasiPeriodicOrbit{D,N}) where {D,N} = N

"""
    jacobi_constant(qpo::QuasiPeriodicOrbit)

Return the average jacobi constant of the states on the invariant curve
"""
function jacobi_constant(qpo::QuasiPeriodicOrbit)
    if typeof(dm(qpo)) <: Cr3bpModel 
        sum(map(x->jacobi_constant(dm(qpo), x0(x)), invariantcurve(qpo)))/length(invariantcurve(qpo)) 
    else
        throw(MethodError(jacobi_constant, qpo))
    end
end

"""
    DGmat(qpo::QuasiPeriodicOrbit)

Return the monodromy matrix of the quasiperiodic orbit
"""
DGmat(qpo::QuasiPeriodicOrbit) = qpo.DG

"""
    eigvals(qpo::QuasiPeriodicOrbit)

Return eigenvalues of the QuasiPeriodicOrbit
"""
LinearAlgebra.eigvals(qpo::QuasiPeriodicOrbit) = copy(qpo.λ)

"""
    eigvecs(qpo::QuasiPeriodicOrbit)

Return eigenvectors of the QuasiPeriodicOrbit
"""
LinearAlgebra.eigvecs(qpo::QuasiPeriodicOrbit) = copy(qpo.V)

"""
    name(qpo::QuasiPeriodicOrbit)

Print the name of the Quasiperiodic Orbit
"""
name(qpo::QuasiPeriodicOrbit) = qpo.name

"""
    family(qpo::QuasiPeriodicOrbit)

Print the family of the Quasiperiodic Orbit
"""
family(qpo::QuasiPeriodicOrbit) = qpo.family

"""
    time2angle(qpo::QuasiPeriodicOrbit)

Convert a ndim time to longitudinal angle
"""
time2angle(qpo::QuasiPeriodicOrbit, T::Real) = 2pi*T/strobetime(qpo)

"""
    angle2time(qpo::QuasiPeriodicOrbit)

Convert a longitudinal angle to ndim time
"""
angle2time(qpo::QuasiPeriodicOrbit, th::Real) = th*strobetime(qpo)/(2pi)

"""
    (qpo::QuasiPeriodicOrbit)(thT; ndtime=false)

Return the invariant curve of the quasiperiodic orbit at time T if ndtime = true, and 
returns the state at thT = 2πT/Period if ndtime = false. thT is the
longitudinal angle on the torus
"""
(qpo::QuasiPeriodicOrbit)(thT; ndtime=false) = ndtime ? invariantcurve(qpo)(thT) : invariantcurve(qpo)(angle2time(qpo, thT))

"""
    (qpo::QuasiPeriodicOrbit)(thT, thrho; ndtime=false)

For the invariant curve `ic` at longitudinal angle `thT` or time `T`, return the 
state on `ic` at latitudinal angle `thrho`
"""
function (qpo::QuasiPeriodicOrbit{dim,N})(thT, thrho; ndtime=false, tol = DEFAULT_CONVERGENCE_TOL) where {dim,N}
    uvec = x2u(qpo(thT; ndtime), xstar(qpo))

    Ut = reshape(uvec, dim, N)

    k = kvec(N)
    Dt = Matrix(transpose(Dmat(N)))


    u = Ut*Dt*Vector(vec(transpose(exp.(im*thrho*k))))

    # Remove imaginary values below a certain tolerance
    maxval, maxind = findmax(x->abs(x), imag(u))
    
    if maxval < tol
        u = real(u)
    else
        println(maxval)
        println(maxind)
        throw(InvalidStateException("Unacceptable imaginary numerical error",:u))
    end
    return u

end
