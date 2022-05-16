# TODO make continuity constraint an abstract type to allow for FullContinuityConstraint
# and PartialContinuityConstraint (or PeriodicityConstraint) concrete types
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                     CONTINUITY CONSTRAINT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
    ContinuityConstraint

Performs x1(T) - x2(T), given x1(0) and x2(T). This constraint
propagates x1 from time 0 to T to perform this operation
"""
struct ContinuityConstraint{D} <: Constraint{D}
    x1::FreeVariable{D,T} where {T}
    x2::FreeVariable{D,T} where {T}
    tof::FreeVariable{1,T} where {T}
    dm::DynamicalModel

    function ContinuityConstraint(x1, x2, T, dm)
        new{length(x1)}(x1,x2,T,dm)
    end
end

"""
    dimension(::ContinuityConstraint{D}) where {D}

Return dimension of the ContinuityConstraint
"""
dimension(::ContinuityConstraint{D}) where {D} = D
Base.length(::ContinuityConstraint{D}) where {D} = D

"""
    x1(cc::ContinuityConstraint)

Return x1
"""
x1(cc::ContinuityConstraint) = cc.x1

"""
    x2(cc::ContinuityConstraint)

Return x2
"""
x2(cc::ContinuityConstraint) = cc.x2

"""
    T(cc::ContinuityConstraint)

Return tof
"""
tof(cc::ContinuityConstraint) = cc.tof

"""
    cctspan(cc::ContinuityConstraint)

Return the tspan to integrate the continuity constraint for
"""
function cctspan(cc::ContinuityConstraint)
    return (0, value(tof(cc))[1])
end

"""
    model(cc::ContinuityConstraint)

Return model
"""
dm(cc::ContinuityConstraint) = cc.dm

"""
    evalconstraint(cc::ContinuityConstraint)

Evaluate the continuity constraint
"""
function evalconstraint(cc::ContinuityConstraint)
    sol = solve(dm(cc), tosvector(x1(cc)), cctspan(cc))
    return sol.u[end]-tovector(x2(cc))
end

"""
    evalconstraint(cc::ContinuityConstraint, X1::FreeVariable, X2::FreeVariable)

Evaluate the continuity constraint
"""
function evalconstraint(cc::ContinuityConstraint, X1::FreeVariable, X2::FreeVariable)
    sol = solve(dm(cc), tosvector(X1), cctspan(cc))
    return sol.u[end]-tovector(X2)
end

"""
    partials(cc::ContinuityConstraint, fv::FreeVariable)

Return the matrix of partial derivatives for the partial of the constraint with
respect to the given free variable
"""
function partials(cc::ContinuityConstraint, fv::FreeVariable)
    if fv == x1(cc)
        # Partial with respect to X1(0)
        return __dCC_dx1
    elseif fv == x2(cc)
        # Partial with respect to X2(0)
        return __dCC_dx2
    elseif  fv == tof(cc) && active(tof(cc))
        # Partial with respect to T
        return __dCC_dt
    else
        # No partial
        return __no_partial
    end
end

"""
    __dCC_dx1

Partial of the continuity constraint with respect to x1(0),
which is the STM phi(0,T)
"""
function __dCC_dx1(cc::ContinuityConstraint)
    sol = tangent_solve(dm(cc), tosvector(x1(cc)), cctspan(cc))
    return sol.u[end][:,2:end]
end

"""
    __dCC_dx2

Partial of the continuity constraint with respect to x2(0),
which is the negative identity matrix
"""
function __dCC_dx2(cc::ContinuityConstraint{D}) where {D}
    return -I(D)
end

"""
    __dCC_dt

Partial of the continuity constraint with respect to T,
which is the derivative of the x1 state at time T
"""
function __dCC_dt(cc::ContinuityConstraint{D}) where {D}
    sol = solve(dm(cc), tosvector(x1(cc)), cctspan(cc))
    return dm(cc)((sol.u[end]))
end

"""
    __no_partial

Partial of the continuity constraint with respect to an
unrelated free variable, which is just a matrix of zeros
"""
function __no_partial(cc::ContinuityConstraint{D}) where {D}
    return zeros(D,D)
end

"""
    Base.show

Overload the show operator to pretty print the ContinuityConstraint to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", cc::ContinuityConstraint)
    print(io, "Continuity Constraint\n")
    print(io, "- Length: $(length(cc))\n")
    print(io, "- X1: $(name(x1(cc)))\n")
    print(io, "- X2: $(name(x2(cc)))\n")
    print(io, "- T: $(name(tof(cc)))\n")
end
