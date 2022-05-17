# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                     CONSTRAINTS
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
using StaticArrays
using LinearAlgebra
"""
    Constraint{D}

Type Parameters
- `D`: dimension of the model
"""
abstract type Constraint{D} end

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
    return zeros(R,C)
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
    tovector(fx::FXVector)

Returns the FXVector as a Vector comprising all the elements of its component FreeVariables
"""
function tovector(fx::FXVector)
    outvec = Vector{Float64}()
    for c in fx
        append!(outvec,evalconstraint(c))
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
