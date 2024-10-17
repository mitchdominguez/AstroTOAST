# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                     CONSTRAINTS
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
using StaticArrays
using LinearAlgebra
using SparseArrays
"""
    Constraint{D}

Type Parameters
- `D`: dimension of the model

Every subtype of Constraint must have a removeinds field
to allow for removing elements from FreeVariables and Constraints
"""
abstract type Constraint{D} end

"""
    dimension(::ContinuityConstraint{D}) where {D}

Return dimension of the ContinuityConstraint
"""
Base.length(::Constraint{D}) where {D} = D

"""
    removeinds(C::Constraint)

Return the indices to remove from the constraint
"""
removeinds(C::Constraint) = copy(C.removeinds)

"""
    evalconstraint(C::Constraint)
"""
evalconstraint(C::Constraint) = throw(MethodError(evalconstraint, C))

"""
    tovector(C::Constraint)

Return the constraint, removing elements
"""
function tovector(C::Constraint)
    return evalconstraint(C)[setdiff(1:end,removeinds(C))]
end

"""
    tofullvector(C::Constraint)

Return the constraint, without removing elements
"""
function tofullvector(C::Constraint)
    return evalconstraint(C)
end



# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                       PARTIAL  
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
"""
    Partial{D}

Abstract type for all partial derivative matrices. The 
type D refers to the number of columns that will be in the 
matrix of partial derivatives
"""
abstract type Partial{D} end

"""
    __NP{D}

Type used when no partial derivative exists.
"""
struct __NP{D} <: Partial{D} end

"""
    (::__NP{C})(cc::Constraint{R})

Function call to output a zero matrix of dimensions RxC
"""
function (::__NP{C})(cc::Constraint{R}) where {R,C}
    return spzeros(full_length(cc),C)
end

# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                       FX Vector
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

struct FXVector{D}
    Cs::Vector{Constraint} 

    # Constructor
    function FXVector(constraints...)
        # Initialize array
        # fv_array = Vector{FreeVariable}(undef,Base.length(constraints))
        fx_array = Vector{Constraint}()
        for i in eachindex(constraints)
            if ~(typeof(constraints[i])<:Constraint)
                throw(MethodError(FXVector, constraints[i]))
            end
            push!(fx_array, constraints[i])
        end
        # new{Base.length(constraints)}(fv_array)
        new{Base.length(fx_array)}(fx_array)
    end
end

""" 
    numels(::XVector)

Returns the number of free variables in the XVector
"""
numels(::FXVector{D}) where {D} = D

"""
    getindex(::FXVector, ...)

Returns the Constraint at index i or in UnitRange r in the given XVector
Note that getindex returns the original FreeVariable objects in the array,
so changes to each fv elsewhere will affect the output here as well
"""
Base.getindex(fx::FXVector, i::Int) = fx.Cs[i]

function Base.getindex(fx::FXVector, r::UnitRange{Int})
    outvec = Vector{Constraint}()
    for i in r
        push!(outvec, fx[i])
    end
    return FXVector(outvec...)
end

"""
    setindex!(fv::FreeVariable{D,T}, ind::Int, newval::T) where {D,T}

Set the value of fv at an index ind to the specified newval
"""
function Base.setindex!(fx::FXVector, newval::T, ind::Int) where {T <: Constraint}
    if 1 <= ind && ind <= length(fx)
        fx.Cs[ind] = newval
    else
        throw(BoundsError(fx,ind)) 
    end
end

"""
    iterate(::FXVector)

Method for iterating through XVectors
"""
Base.iterate(fx::FXVector, state=1) = state>numels(fx) ? nothing : (fx[state], state+1)

"""
    length(fx::FXVector)

Returns the length of the FXVector if it is output as a Vector
"""
function Base.length(fx::FXVector)
    els = 0
    for c in fx
        els = els + length(c)
    end
    return els
end

"""
    removeinds(fx::FXVector)

Output a vector containing the indices of removed elements in the full XVector
"""
function removeinds(fx::FXVector)
    ind = 0
    outvec = Vector{Int}()
    for i = 1:numels(fx)
        append!(outvec, removeinds(fx[i]).+ind)
        ind += full_length(fx[i])
    end
    return outvec
end

"""
    tovector(fx::FXVector)

Returns the FXVector as a Vector comprising all the elements of its component FreeVariables.
This function removes the elements specified in removeinds
"""
function tovector(fx::FXVector)
    outvec = Vector{Float64}()
    for c in fx
        append!(outvec,tovector(c))
    end
    return outvec
end

"""
    tofullvector(fx::FXVector)

Returns the FXVector as a Vector comprising all the elements of its component FreeVariables
"""
function tofullvector(fx::FXVector)
    outvec = Vector{Float64}()
    for c in fx
        append!(outvec,tofullvector(c))
    end
    return outvec
end

"""
    tovector(fx::FXVector, xvorig::XVector, xveval::XVector)

Returns the FXVector as a Vector comprising all the elements of its component FreeVariables
This method of `tovector` assumes that the constraints within fx are all defined with
respect to free variables contained within xvorig. To evaluate constraints,
this method will update! xvorig to the values contained in xveval, evaluate
the constraint, and then update! xvorig back to its original state. 
"""
function tovector(fx::FXVector, xvorig::XVector, xveval::Vector{FreeVariable})
    xvorig_copy = copy(xvorig)
    update!(xvorig, xveval)
    outvec = tovector(fx)
    update!(xvorig, xvorig_copy)
    return outvec
end

function tovector(fx::FXVector, xvorig::XVector, xveval::XVector)
    tovector(fx, xvorig, xveval.FVs)
end

"""
    tosvector(fx::FXVector)

Returns the FXVector as an SVector comprising all the elements of its component FreeVariables
"""
function tosvector(fx::FXVector)
    SVector(Tuple(tovector(fx)))
end

"""
    norm(fx::FXVector)

Calculate the L2 norm of the constraint vector
"""
LinearAlgebra.norm(fx::FXVector) = LinearAlgebra.norm(tovector(fx))
LinearAlgebra.norm(fx::FXVector, xvorig::XVector, xveval::XVector) = LinearAlgebra.norm(tovector(fx, xvorig, xveval))

"""
    Base.show

Overload the show operator to pretty print the FXVector to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", fx::FXVector)
    print(io, "FXVector:\n")
    print(io, "- Number of Constraints: $(numels(fx))\n")
    print(io, "- Length: $(length(fx))\n")
    print(io, "---\n")
    print(io, "Constraints:\n")
    for c in fx
        # print(io, "- Name: $(name(fv)), Length: $(length(fv))\n")
        show(io, "text/plain", c)
    end
end
