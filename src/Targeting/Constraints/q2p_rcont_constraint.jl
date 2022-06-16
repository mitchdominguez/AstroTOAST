# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                       QPO-PO POSITION CONTINUITY CONSTRAINT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
    Q2P_PositionContinuityConstraint

Invariance constraint for a quasiperiodic orbit on a 2-dimensional torus
"""
struct Q2P_PositionContinuityConstraint{D} <: Constraint{D}
    qpo::QuasiPeriodicOrbit{D1,N} where {D1,N}
    po::PeriodicOrbit{D1}  where {D1}
    thT::FreeVariable
    thR::FreeVariable
    tau::FreeVariable
    removeinds::Vector{Int}

    function Q2P_PositionContinuityConstraint(qpo, po, thT, thR, tau)
        if dimension(qpo) != dimension(po)
            throw(ErrorException("Dimension of PO $(dimension(po))and QPO $(dimension(qpo)) are incompatible"))
        end

        if dm(qpo) != dm(po)
            throw(ErrorException("DynamicalModel of PO $(dm(po))and QPO $(dm(qpo)) are incompatible"))
        end

        if !(full_length(thT) == full_length(thR) == full_length(tau) == 1)
            throw(ErrorException("thT, thR, and tau must all have full_length 1"))
        end

        new{3}(qpo, po, thT, thR, tau, [])
    end
end


#################################################################
# Basic functions related to the Q2P_PositionContinuityConstraint
#################################################################

"""
    dm(q2p::Q2P_PositionContinuityConstraint)

Return the dynamical model of the InvarianceConstraint
"""
dm(q2p::Q2P_PositionContinuityConstraint) = dm(q2p_po(q2p))

"""
    full_length(q2p::Q2P_PositionContinuityConstraint)

Return the full_length of the Q2P_PCC
"""
full_length(q2p::Q2P_PositionContinuityConstraint) = length(q2p)

"""
    q2p_qpo(q2p::Q2P_PositionContinuityConstraint)

Return the QPO part of the constraint
"""
q2p_qpo(q2p::Q2P_PositionContinuityConstraint) = q2p.qpo

"""
    q2p_po(q2p::Q2P_PositionContinuityConstraint)

Return the PO part of the constraint
"""
q2p_po(q2p::Q2P_PositionContinuityConstraint) = q2p.po

"""
    q2p_tht(q2p::Q2P_PositionContinuityConstraint)

Return the FreeVariable corresponding to longitudinal angle of the QPO
"""
q2p_tht(q2p::Q2P_PositionContinuityConstraint) = q2p.thT

"""
    q2p_thr(q2p::Q2P_PositionContinuityConstraint)

Return the FreeVariable corresponding to latitudinal angle of the QPO
"""
q2p_thr(q2p::Q2P_PositionContinuityConstraint) = q2p.thR

"""
    q2p_tau(q2p::Q2P_PositionContinuityConstraint)

Return the FreeVariable corresponding to longitudinal angle of the PO
"""
q2p_tau(q2p::Q2P_PositionContinuityConstraint) = q2p.tau



#################################################################
######### Solving the Q2P_PositionContinuityConstraint ##########
#################################################################

"""
    evalconstraint(q2p::Q2P_PositionContinuityConstraint)

Evaluate the constraint
"""
function evalconstraint(q2p::Q2P_PositionContinuityConstraint)
    po = q2p_po(q2p)
    qpo = q2p_qpo(q2p)
    thT = tofullvector(q2p_tht(q2p))[1]
    thR = tofullvector(q2p_thr(q2p))[1]
    tau = tofullvector(q2p_tau(q2p))[1]

    return (qpo(thT,thR)-po(tau))[1:3]
end

"""
    partials(q2p::Q2P_PositionContinuityConstraint, fv::FreeVariable)

Return the matrix of partial derivatives for the partial of the constraint with
respect to the given free variable
"""
function partials(q2p::Q2P_PositionContinuityConstraint, fv::FreeVariable{D,T}) where {D,T}
    if fv == q2p_tht(q2p) && active(q2p_tht(q2p))
        # Partial with respect to the invariant curve states
        return __dQ2P_dthT{D}()

    elseif fv == q2p_thr(q2p) && active(q2p_thr(q2p))
        # Partial with respect to stroboscopic time
        return __dQ2P_dthR{D}()

    elseif fv == q2p_tau(q2p) && active(q2p_tau(q2p))
        # Partial with respect to rotation angle
        return __dQ2P_dtau{D}()

    else
        # No partial
        return __NP{full_length(fv)}()
    end
end

"""
    __dQ2P_dthT

Partial of the Q2P_PositionContinuityConstraint with respect to the longitudinal angle
"""
struct __dQ2P_dthT{D} <: Partial{D} end
function (::__dQ2P_dthT{C})(q2p::Q2P_PositionContinuityConstraint{R}) where {R,C}
    # Define q2p variables
    qpo = q2p_qpo(q2p)
    model = dm(q2p)
    thT = tofullvector(q2p_tht(q2p))[1]
    thR = tofullvector(q2p_thr(q2p))[1]

    # States
    x_qpo = qpo(thT,thR)

    # Angular frequencies
    ωT = 2π/strobetime(qpo)
    ωR = rotationangle(qpo)/strobetime(qpo)

    # Calculate partials
    dudthR = __dQ2P_dthR{R}()(q2p)
    dudthT = (1/ωT)*(model(x_qpo)[1:3] - ωR*dudthR)

    return dudthT[1:3]
end

"""
    __dQ2P_dthR

Partial of the Q2P_PositionContinuityConstraint with respect to the latitudinal angle
"""
struct __dQ2P_dthR{D} <: Partial{D} end
function (::__dQ2P_dthR{C})(q2p::Q2P_PositionContinuityConstraint{R}) where {R,C}
    # Define q2p variables
    po = q2p_po(q2p)
    qpo = q2p_qpo(q2p)
    thT = tofullvector(q2p_tht(q2p))[1]
    thR = tofullvector(q2p_thr(q2p))[1]
    model = dm(q2p)
    N = numels(qpo)
    k = kvec(N)
    kdiag = Diagonal(vec(k))
    D = Dmat(N)
    U = transpose(reshape(x2u(qpo(thT), xstar(qpo,thT)), dimension(qpo), N))

    # Partials
    dudthR = ignore_imag(vec(im*exp.(im*thR*k)*kdiag*D*U), DEFAULT_CONVERGENCE_TOL)

    return dudthR[1:3]
end

"""
    __dQ2P_dtau

Partial of the Q2P_PositionContinuityConstraint with respect to the tau angle
"""
struct __dQ2P_dtau{D} <: Partial{D} end
function (::__dQ2P_dtau{C})(q2p::Q2P_PositionContinuityConstraint{R}) where {R,C}
    # Define q2p variables
    po = q2p_po(q2p)
    tau = tofullvector(q2p_tau(q2p))[1]
    model = dm(q2p)

    return -model(po(tau))[1:3]
end

"""
    Base.show

Overload the show operator to pretty print the PeriodicOrbit to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", q2p::Q2P_PositionContinuityConstraint)
    print(io, "Q2P_PositionContinuityConstraint:\n")
    print(io, "\n===================================\n")
    show(io, "text/plain", q2p_qpo(q2p))
    print(io, "\n===================================\n")
    show(io, "text/plain", q2p_po(q2p))
    print(io, "\n- FreeVariables:\n===================================\n")
    print(io, "thT = $(tofullvector(q2p_tht(q2p))[1])\n")
    print(io, "thR = $(tofullvector(q2p_thr(q2p))[1])\n")
    print(io, "tau = $(tofullvector(q2p_tau(q2p))[1])\n")
end
