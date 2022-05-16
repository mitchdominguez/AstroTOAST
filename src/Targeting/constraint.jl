# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                     CONSTRAINTS
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
using StaticArrays
"""
    Constraint{D}

Type Parameters
- `D`: dimension of the model
"""
abstract type Constraint{D} end

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
    tosvector(fx::FXVector)

Returns the FXVector as an SVector comprising all the elements of its component FreeVariables
"""
function tosvector(fx::FXVector)
    SVector(Tuple(tovector(fx)))
end

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
