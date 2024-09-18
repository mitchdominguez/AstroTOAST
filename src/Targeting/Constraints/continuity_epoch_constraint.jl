# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                     CONTINUITY CONSTRAINT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
# TODO should I make a new field called solvefunction which I can use to call different ODE solvers, including MATLAB?
# whatever this is, it would still need the same three argument interface as what currently exists
"""
    ContinuityEpochConstraint

Performs x1(T) - x2(T), given x1(0) and x2(T). This constraint
propagates x1 from time 0 to T to perform this operation
"""
struct ContinuityEpochConstraint{D} <: Constraint{D}
    x1::FreeVariable{D1,T} where {D1,T} # TODO D for X1 and X2 should be allowed to differ
    x2::FreeVariable{D2,T} where {D2,T}
    epoch1::FreeVariable{1,T} where {T}
    epoch2::FreeVariable{1,T} where {T}
    tof::FreeVariable{N,T} where {N,T}
    dm::DynamicalModel
    removeinds::Vector{Int}
    solvefunction::Function
    tangentsolvefunction::Function

    # function ContinuityEpochConstraint(x1, x2, T, dm, rminds = Vector{Int}()::Union{Int, AbstractVector{Int}})
    function ContinuityEpochConstraint(x1, x2, T, dm; includeinds = Vector(1:dimension(dm)), solvefunction=solve, tangentsolvefunction=tangent_solve)


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
    solvefunction(cc::ContinuityEpochConstraint)

Return the function that is used to propagate the states within the constraint
"""
solvefunction(cc::ContinuityEpochConstraint) = cc.solvefunction

"""
    tangentsolvefunction(cc::ContinuityEpochConstraint)

Return the function that is used to propagate the states within the constraint
"""
tangentsolvefunction(cc::ContinuityEpochConstraint) = cc.tangentsolvefunction

"""
    full_length(cc::ContinuityEpochConstraint)

Return the full length of the ContinuityEpochConstraint, without
removing elements
"""
full_length(cc::ContinuityEpochConstraint) = full_length(x1(cc))

"""
    x1(cc::ContinuityEpochConstraint)

Return x1
"""
x1(cc::ContinuityEpochConstraint) = cc.x1

"""
    x2(cc::ContinuityEpochConstraint)

Return x2
"""
x2(cc::ContinuityEpochConstraint) = cc.x2

"""
    T(cc::ContinuityEpochConstraint)

Return tof
"""
tof(cc::ContinuityEpochConstraint) = cc.tof

"""
    get_epoch1(cc::ContinuityEpochConstraint)

Return get_epoch1
"""
get_epoch1(cc::ContinuityEpochConstraint) = cc.epoch1

"""
    get_epoch2(cc::ContinuityEpochConstraint)

Return get_epoch2
"""
get_epoch2(cc::ContinuityEpochConstraint) = cc.epoch2

"""
    cctspan(cc::ContinuityEpochConstraint)

Return the tspan to integrate the continuity constraint for
"""
function cctspan(cc::ContinuityEpochConstraint)
    e1 = get_epoch1(cc)
    return (e1, e1+fullvalue(tof(cc))[1])
end

"""
    model(cc::ContinuityEpochConstraint)

Return model
"""
dm(cc::ContinuityEpochConstraint) = cc.dm

"""
    evalconstraint(cc::ContinuityEpochConstraint)

Evaluate the continuity constraint
"""
function evalconstraint(cc::ContinuityEpochConstraint)
    F16 = solvefunction(cc)(dm(cc), tofullsvector(x1(cc)), cctspan(cc)).u[end]-tofullvector(x2(cc))
    F7 = get_epoch1(cc) + tof(cc) - get_epoch2(cc)

    return vcat(F16...,F7)
end

# """
    # evalconstraint(cc::ContinuityEpochConstraint, X1::FreeVariable, X2::FreeVariable)

# Evaluate the continuity constraint
# """
# function evalconstraint(cc::ContinuityEpochConstraint, X1::FreeVariable, X2::FreeVariable)
    # sol = solvefunction(cc)(dm(cc), tofullsvector(X1), cctspan(cc))
    # return sol.u[end]-tofullvector(X2)
# end

"""
    partials(cc::ContinuityEpochConstraint, fv::FreeVariable)

Return the matrix of partial derivatives for the partial of the constraint with
respect to the given free variable
"""
function partials(cc::ContinuityEpochConstraint, fv::FreeVariable{D,T}) where {D,T}
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

    elseif fv == get_epoch1(cc) && active(get_epoch1(cc))
        # Partial with respect to Epoch 1 (epoch of first patch point)
        return __dCC_depoch1{D}()

    elseif fv == get_epoch2(cc) && active(get_epoch2(cc))
        # Partial with respect to Epoch 2 (epoch of first patch point)
        return __dCC_depoch2{D}()

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
function (::__dCC_dx_ss{C})(cc::ContinuityEpochConstraint{R}) where {R,C}
    sol = tangentsolvefunction(cc)(dm(cc), tofullsvector(x1(cc)), cctspan(cc))

    # println("TSPAN")
    # println(cctspan(cc))
    # println("-----")
    # show(stdout,"text/plain",sol.u[end][:,2:end])
    # println("")
    return sol.u[end][:,2:(end-1)] - I(full_length(x1(cc)))
end

"""
    __dCC_dx1

Partial of the continuity constraint with respect to x1(0),
which is the STM phi(0,T)
"""
struct __dCC_dx1{D} <: Partial{D} end
function (::__dCC_dx1{C})(cc::ContinuityEpochConstraint{R}) where {R,C}
    sol = tangentsolvefunction(cc)(dm(cc), tofullsvector(x1(cc)), cctspan(cc))

    # println("TSPAN")
    # println(cctspan(cc))
    # println("-----")
    # show(stdout,"text/plain",sol.u[end][:,2:end])
    # println("")
    return vcat(sol.u[end][:,2:(end-1)], zeros(1,dimension(dm(cc))))
end

"""
    __dCC_dx2

Partial of the continuity constraint with respect to x2(0),
which is the negative identity matrix
"""
struct __dCC_dx2{D} <: Partial{D} end
function (::__dCC_dx2{C})(cc::ContinuityEpochConstraint{R}) where {R,C}
    return vcat(-I(C), zeros(1,dimension(dm(cc))))
end

"""
    __dCC_dt

Partial of the continuity constraint with respect to T,
which is the derivative of the x1 state at time T
"""
struct __dCC_dt{D} <: Partial{D} end
function (::__dCC_dt{C})(cc::ContinuityEpochConstraint{R}) where {R,C}
    sol = solvefunction(cc)(dm(cc), tofullsvector(x1(cc)), cctspan(cc))
    return vcat(dm(cc)((sol.u[end]), tof(cc)[1] + get_epoch1(cc)[1]), 1.0)
end

"""
    __dCC_depoch1

Partial of the continuity constraint with respect to the epoch of the first patch point
"""
struct __dCC_depoch1{D} <: Partial{D} end
function (::__dCC_depoch1{C})(cc::ContinuityEpochConstraint{R}) where {R,C}
    epochpartial = tangentsolvefunction(cc)(dm(cc), tofullsvector(x1(cc)), cctspan(cc)).u[end][:,end]
    return vcat(epochpartial, 1.0)
end

"""
    __dCC_depoch2

Partial of the continuity constraint with respect to the epoch of the first patch point
"""
struct __dCC_depoch2{D} <: Partial{D} end
function (::__dCC_depoch2{C})(cc::ContinuityEpochConstraint{R}) where {R,C}
    return vcat(zeros(dimension(dm(cc))), -1.0) # TODO check that this is the right dimension to be returning
end

"""
    Base.show

Overload the show operator to pretty print the ContinuityEpochConstraint to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", cc::ContinuityEpochConstraint)
    print(io, "Continuity Constraint\n")
    print(io, "- Length: $(length(cc))\n")
    print(io, "- X1: $(name(x1(cc)))\n")
    print(io, "- X2: $(name(x2(cc)))\n")
    print(io, "- T: $(name(tof(cc)))\n")
end
