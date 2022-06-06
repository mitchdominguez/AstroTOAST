# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                  INVARIANCE CONSTRAINT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
    InvarianceConstraint2D

Invariance constraint for a quasiperiodic orbit on a 2-dimensional torus
"""
struct InvarianceConstraint2D{D} <: Constraint{D}
    u0::FreeVariable{D1,T} where {D1,T}
    xstar::Vector{Float64}
    T::FreeVariable{D3,T} where {D3,T}
    ρ::FreeVariable{D4,T} where {D4,T}
    dm::DynamicalModel
    removeinds::Vector{Int}

    function InvarianceConstraint2D(u0, xstar, T, ρ, dm, rminds = Vector{Int}()::Union{Int, AbstractVector{Int}})
        dim = dimension(dm)

        if full_length(u0)%dim == 0 && length(xstar) == dim
            N = Int(full_length(u0)/dim)

            # Check that N is odd
            if iseven(N)
                throw(ExceptionError("There must be an odd number of nodes in the invariant curve"))
            end

            # Check that rho and T are both 1 dimensional
            if full_length(ρ) != 1 || full_length(T) != 1
                throw(DimensionMismatch("ρ and T must have full length 1"))
            end

            rmlength = Base.length(rminds)
            vallength = full_length(u0)

            # Check for correct bounds on removeinds
            if !all(rminds.<=vallength) || !all(rminds.>=1) || rmlength>vallength
                throw(BoundsError(value, removeinds))
            end
            
            # Calculate dimension of continuity constraint
            rmvec = Vector{Int}()
            append!(rmvec, rminds)
            
            icdimension = vallength-length(rmvec) # D in FreeVariable{D,T}

            new{icdimension}(u0, xstar, T, ρ, dm, rmvec)
        else
            throw(DimensionMismatch("dimension of u0, xstar, and model are incompatible"))
        end
    end
end

#########################################################
# Basic functions related to the InvarianceConstraint2d
#########################################################

"""
    u0(ic::InvarianceConstraint2D)

Return the invariant curve FreeVariable
"""
u0(ic::InvarianceConstraint2D) = ic.u0

"""
    u0vec(ic::InvarianceConstraint2D)

Return the invariant curve points as a vecotr
"""
u0vec(ic::InvarianceConstraint2D) = tofullvector(u0(ic))

"""
    xstar(ic::InvarianceConstraint2D)

Return the fixed point of the invariant curve
"""
xstar(ic::InvarianceConstraint2D) = ic.xstar

"""
    strobetime(ic::InvarianceConstraint2D)

Return the stroboscopic period of the torus to be targeted
"""
strobetime(ic::InvarianceConstraint2D) = tofullvector(ic.T)[1]

"""
    rotationangle(ic::InvarianceConstraint2D)

Return the angle by which the invariant curve will rotate after one period T
"""
rotationangle(ic::InvarianceConstraint2D) = tofullvector(ic.ρ)[1]

"""
    full_length(ic::InvarianceConstraint2D)

Return the full length of the ContinuityConstraint, without
removing elements
"""
full_length(ic::InvarianceConstraint2D) = full_length(u0(ic))

"""
    dm(ic::InvarianceConstraint2D)

Return the dynamical model of the InvarianceConstraint
"""
dm(ic::InvarianceConstraint2D) = ic.dm

############################################################
# Utility functions for solving the invariance constraint
############################################################

"""
    N(ic::InvarianceConstraint2D)

Return the number of nodes in the invariance constraint
"""
numels(ic::InvarianceConstraint2D) = Int(full_length(ic)/dimension(dm(ic)))

"""
    thvec(ic::InvarianceConstraint2D)

Generate the vector of angles (in radians) used in the DFT
"""
thvec(N::Int) = 2pi/N*(0:N-1)'

"""
    k(ic::InvarianceConstraint2D)

Generate k vector used in the DFT
"""
kvec(N::Int) = Vector{Int}(-(N-1)/2:(N-1)/2)'

"""
    Dmat(N::Int)

Calculate D matrix used in the DFT
"""
Dmat(N::Int) = 1/N*exp.(-im*kvec(N)'*thvec(N))

"""
    Qmat(N::Int, ρ::Float64)

Calculate Q matrix used in the DFT
"""
Qmat(N::Int, ρ::Float64) = Diagonal(vec(exp.(im*kvec(N)*ρ)))

"""
    invariant_rotation(ic::InvarianceConstraint2D, ρ::Float64)

Calculate the rotation matrix for rotating an angle `ρ` about the
invariant curve
"""
function invariant_rotation(ic::InvarianceConstraint2D, ρ::Float64; tol=1e-15)
    N = numels(ic)
    D = Dmat(N)
    Q = Qmat(N, ρ)
    R_nr = D\Q*D

    # Remove imaginary values below a certain tolerance
    maxval, maxind = findmax(x->abs(x), imag(R_nr))
    
    if maxval < tol
        R_nr = real(R_nr)
    else
        throw(ErrorException("Unacceptable imaginary numerical error"))
    end
end

"""
    propagate_invariant_curve(ic::InvarianceConstraint2D)

Propagate each node on the invariant curve for one stroboscopic
time period
"""
function propagate_invariant_curve(ic::InvarianceConstraint2D)
    N = numels(ic)
    u0 = u0vec(ic)
    utr = similar(u0)

    for i = 1:N
        q0 = u0[6*i-5:6*i] + xstar(ic) # Full state including fixed point
        sol = solve(dm(ic), q0, strobetime(ic)) # Propagate state forwards
        utr[6*i-5:6*i] = sol[end]-xstar(ic) # Remove fixed point from propagated answer
    end

    return utr
end

"""
    evalconstraint(ic::InvarianceConstraint2D)

Evaluate the invariance constraint
"""
function evalconstraint(ic::InvarianceConstraint2D)
    # # Calculate dQdrho
    # dQdrho = Diagonal(vec(k))*(-im*Q)
    # # Calculate partial of rotation matrix wrt rho
    # dRdrho = D\dQdrho*D

    # Rotation matrix
    R_nr = invariant_rotation(ic, -rotationangle(ic))
    
    utr = kron(R_nr, I(6))*propagate_invariant_curve(ic)
    
    return utr - tofullvector(u0(ic))
end
