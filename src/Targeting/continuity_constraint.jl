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
    tof::FreeVariable{N,T} where {N,T}
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
    return (0, fullvalue(tof(cc))[1])
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
    sol = solve(dm(cc), tofullsvector(x1(cc)), cctspan(cc))
    return sol.u[end]-tovector(x2(cc)) ### TODO this needs to be changed
end

"""
    evalconstraint(cc::ContinuityConstraint, X1::FreeVariable, X2::FreeVariable)

Evaluate the continuity constraint
"""
function evalconstraint(cc::ContinuityConstraint, X1::FreeVariable, X2::FreeVariable)
    sol = solve(dm(cc), tofullsvector(X1), cctspan(cc))
    # deleteat!(copy(fv.value),removeinds(fv))
    return sol.u[end]-tovector(X2) ### TODO this needs to be changed
end

"""
    partials(cc::ContinuityConstraint, fv::FreeVariable)

Return the matrix of partial derivatives for the partial of the constraint with
respect to the given free variable
"""
function partials(cc::ContinuityConstraint, fv::FreeVariable{D,T}) where {D,T}
    if fv == x1(cc)
        # Partial with respect to X1(0)
        return __dCC_dx1{D}() # TODO make this have a constructor that can take in fv, cc
    elseif fv == x2(cc)
        # Partial with respect to X2(0)
        return __dCC_dx2{D}()
    elseif  fv == tof(cc) && active(tof(cc))
        # Partial with respect to T
        return __dCC_dt{D}()
    else
        # No partial
        return __NP{D}()
    end
end

"""
    __dCC_dx1

Partial of the continuity constraint with respect to x1(0),
which is the STM phi(0,T)
"""
struct __dCC_dx1{D} <: Partial{D} end
function (::__dCC_dx1{C})(cc::ContinuityConstraint{R}) where {R,C}
    sol = tangent_solve(dm(cc), tofullsvector(x1(cc)), cctspan(cc))
    return sol.u[end][:,2:end]
end

"""
    __dCC_dx2

Partial of the continuity constraint with respect to x2(0),
which is the negative identity matrix
"""
struct __dCC_dx2{D} <: Partial{D} end
function (::__dCC_dx2{C})(cc::ContinuityConstraint{R}) where {R,C}
    return -I(C)
end

"""
    __dCC_dt

Partial of the continuity constraint with respect to T,
which is the derivative of the x1 state at time T
"""
struct __dCC_dt{D} <: Partial{D} end
function (::__dCC_dt{C})(cc::ContinuityConstraint{R}) where {R,C}
    sol = solve(dm(cc), tofullsvector(x1(cc)), cctspan(cc))
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
