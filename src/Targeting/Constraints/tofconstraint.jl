# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                TIME OF FLIGHT CONSTRAINT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
    TOFConstraint

Applies the constraint ΣT_i - T_d = 0, where T_i are FreeVariables for times of flight,
and T_d is the desired total time of flight
"""
struct TOFConstraint{D} <: Constraint{D}
    T::Vector{FreeVariable}
    Td::Real
    removeinds::Vector{Int}

    # Constructor
    function TOFConstraint(Td, T...)

        # Check that Td is nonnegative
        if Td < 0
            throw(InvalidStateException("Td must be nonnegative", :Td))
        end

        # Check that all T are valid FreeVariables
        Tvec = Vector{FreeVariable}()
        for Ti in T
            # Check that all elements of T are FreeVariables
            if typeof(Ti)<:FreeVariable

                # Check that all FreeVariables in T are of full_length 1
                if full_length(Ti) == 1
                    push!(Tvec, Ti)
                else
                    throw(DimensionMismatch("All FreeVariables passed into the
                                            TOFConstraint constructor must be
                                            of full_length = 1"))
                end # Check Ti length
            else
                throw(MethodError(TOFConstraint, Ti))
            end # Check freevariable

        end # for Ti in T

        new{1}(Tvec, Td, Vector{Int}([]))
    end

    TOFConstraint() = new{1}()
end

"""
    td(tofc::TOFConstraint)

Return the desired time of flight
"""
td(tofc::TOFConstraint) = tofc.Td

"""
    tvec(tofc::TOFConstraint)

Return the Vector of TOF FreeVariables
"""
tvec(tofc::TOFConstraint) = tofc.T

"""
    full_length(tofc::TOFConstraint)

Return the full length of the TOFConstraint
"""
full_length(tofc::TOFConstraint) = length(tofc)

"""
    evalconstraint(tofc::TOFConstraint)

Evaluate the continuity constraint
"""
function evalconstraint(tofc::TOFConstraint)
    T_tot = 0
    T = tvec(tofc)
    for Ti in T
        T_tot+=tofullvector(Ti)[1]
    end
    return [T_tot - td(tofc)]
end

"""
    evalconstraint(cc::TOFConstraint, Td::Real, T::XVector)

Evaluate the continuity constraint
"""
function evalconstraint(cc::TOFConstraint, Td::Real, T::XVector)
    temptofc = TOFConstraint(copy(Td), copy(T)...)
    return [evalconstraint(temptofc)]
end

"""
    partials(tofc::TOFConstraint, fv::FreeVariable)

Return the matrix of partial derivatives for the partial of the constraint with
respect to the given free variable
"""
function partials(tofc::TOFConstraint, fv::FreeVariable{D,T}) where {D,T}
    if fv in tvec(tofc)
        return __dTOFC_dTi{full_length(fv)}()
    else
        # No partial
        # return __NP{D}()
        return __NP{full_length(fv)}()
    end
end

"""
    __dTOFC_dTi

Partial of the time of flight constraint with respect to Ti, which
is a FreeVariable included in the TOF constraint
"""
struct __dTOFC_dTi{D} <: Partial{D} end
function (::__dTOFC_dTi{C})(cc::TOFConstraint{R}) where {R,C}
    return ones(R,C)
end
