# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                     CONTINUITY CONSTRAINT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
# TODO should I make a new field called solvefunction which I can use to call different ODE solvers, including MATLAB?
# whatever this is, it would still need the same three argument interface as what currently exists
"""
    ContinuityConstraint

Performs x1(T) - x2(T), given x1(0) and x2(T). This constraint
propagates x1 from time 0 to T to perform this operation
"""
struct ContinuityConstraint{D} <: Constraint{D}
    x1::FreeVariable{D1,T} where {D1,T} # TODO D for X1 and X2 should be allowed to differ
    x2::FreeVariable{D2,T} where {D2,T}
    tof::FreeVariable{N,T} where {N,T}
    dm::DynamicalModel
    removeinds::Vector{Int}
    solvefunction::Function
    tangentsolvefunction::Function

    # function ContinuityConstraint(x1, x2, T, dm, rminds = Vector{Int}()::Union{Int, AbstractVector{Int}})
    function ContinuityConstraint(x1, x2, T, dm; includeinds = Vector(1:dimension(dm)), solvefunction=solve, tangentsolvefunction=tangent_solve)


        if full_length(x1) == full_length(x2) == dimension(dm)

            # Find which indices to remove. Note that because of how
            # setdiff works, it is impossible to make rminds have 
            # indices that are outside the bounds of x1 or x2
            rminds = setdiff(Vector(1:dimension(dm)), includeinds)

            rmlength = Base.length(rminds)
            vallength = full_length(x1)

            # Check for correct bounds on removeinds
            if !all(rminds.<=vallength) || !all(rminds.>=1) || rmlength>vallength
                throw(BoundsError(value, removeinds))
            end
            
            # Calculate dimension of continuity constraint
            #   Combine removeinds(X1), removeinds(X2), and removeinds(cc)
            # rmvec = sort(union(removeinds(x1),removeinds(x2),rminds))
            rmvec = Vector{Int}()
            append!(rmvec, rminds)
            
            ccdimension = vallength-length(rmvec) # D in FreeVariable{D,T}

            new{ccdimension}(x1,x2,T,dm,rmvec,solvefunction,tangentsolvefunction)
        else
            throw(DimensionMismatch("X1 and X2 must have the same full length"))
        end
    end
end

"""
    solvefunction(cc::ContinuityConstraint)

Return the function that is used to propagate the states within the constraint
"""
solvefunction(cc::ContinuityConstraint) = cc.solvefunction

"""
    tangentsolvefunction(cc::ContinuityConstraint)

Return the function that is used to propagate the states within the constraint
"""
tangentsolvefunction(cc::ContinuityConstraint) = cc.tangentsolvefunction

"""
    full_length(cc::ContinuityConstraint)

Return the full length of the ContinuityConstraint, without
removing elements
"""
full_length(cc::ContinuityConstraint) = full_length(x1(cc))

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
    sol = solvefunction(cc)(dm(cc), tofullsvector(x1(cc)), cctspan(cc))
    return sol.u[end]-tofullvector(x2(cc))
end

"""
    evalconstraint(cc::ContinuityConstraint, X1::FreeVariable, X2::FreeVariable)

Evaluate the continuity constraint
"""
function evalconstraint(cc::ContinuityConstraint, X1::FreeVariable, X2::FreeVariable)
    sol = solvefunction(cc)(dm(cc), tofullsvector(X1), cctspan(cc))
    return sol.u[end]-tofullvector(X2)
end

"""
    partials(cc::ContinuityConstraint, fv::FreeVariable)

Return the matrix of partial derivatives for the partial of the constraint with
respect to the given free variable
"""
function partials(cc::ContinuityConstraint, fv::FreeVariable{D,T}) where {D,T}
    if fv == x1(cc) == x2(cc)
        # Constraint takes the form X1(T) - X1(0)
        # Note that this is the single shooter case
        return __dCC_dx_ss{D}()
    elseif fv == x1(cc) 
        # Partial with respect to X1(T)
        return __dCC_dx1{D}() 
    elseif fv == x2(cc)
        # Partial with respect to X2(0)
        return __dCC_dx2{full_length(fv)}()
    elseif  fv == tof(cc) && active(tof(cc))
        # Partial with respect to T
        return __dCC_dt{D}()
    else
        # No partial
        # return __NP{D}()
        return __NP{full_length(fv)}()
    end
end

"""
    __dCC_dx_ss

Partial of the continuity constraint with respect to x1(0),
which is the STM phi(0,T)
"""
struct __dCC_dx_ss{D} <: Partial{D} end
function (::__dCC_dx_ss{C})(cc::ContinuityConstraint{R}) where {R,C}
    sol = tangentsolvefunction(cc)(dm(cc), tofullsvector(x1(cc)), cctspan(cc))

    # println("TSPAN")
    # println(cctspan(cc))
    # println("-----")
    # show(stdout,"text/plain",sol.u[end][:,2:end])
    # println("")
    return sol.u[end][:,2:end] - I(full_length(x1(cc)))
end

"""
    __dCC_dx1

Partial of the continuity constraint with respect to x1(0),
which is the STM phi(0,T)
"""
struct __dCC_dx1{D} <: Partial{D} end
function (::__dCC_dx1{C})(cc::ContinuityConstraint{R}) where {R,C}
    sol = tangentsolvefunction(cc)(dm(cc), tofullsvector(x1(cc)), cctspan(cc))

    # println("TSPAN")
    # println(cctspan(cc))
    # println("-----")
    # show(stdout,"text/plain",sol.u[end][:,2:end])
    # println("")
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
    sol = solvefunction(cc)(dm(cc), tofullsvector(x1(cc)), cctspan(cc))
    return dm(cc)((sol.u[end]))
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
