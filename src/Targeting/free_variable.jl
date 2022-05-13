# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                    FREE VARIABLE
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
import Base.getindex
import Base.length
import Base.iterate
"""
    FreeVariable{D}

"""
struct FreeVariable{D,T}
    name::String
    value::T
    
    # Constructor
    function FreeVariable(name, value)
        new{Base.length(value), typeof(value)}(name, value) 
    end
end

"""
    dimension(::FreeVariable)

Return the dimension of the free variable, i.e. the number of state vector elements
"""
Base.length(::FreeVariable{D}) where {D} = D

"""
    name(fv::FreeVariable)

Return the string name of the free variable
"""
name(fv::FreeVariable) = fv.name

"""
    value(fv::FreeVariable)

Return the value of the free variable
"""
value(fv::FreeVariable) = fv.value

"""
    getindex(::FreeVariable, i::Int)

Returns the element of the FreeVariable at index i in the given XVector
"""
Base.getindex(fv::FreeVariable, i::Int) = value(fv)[i]

"""
    iterate(::FreeVariable)

Method for iterating through FreeVariables
"""
Base.iterate(fv::FreeVariable, state=1) = state>length(fv) ? nothing : (fv[state], state+1)

# TODO set?


# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                       XVector
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
"""
    XVector
"""
struct XVector{D}
    FVs::Vector{FreeVariable} 
    
    # Constructor
    function XVector(free_variables...)
        # Initialize array
        fv_array = Vector{FreeVariable}(undef,Base.length(free_variables))
        for i in eachindex(free_variables)
            if ~(typeof(free_variables[i])<:FreeVariable)
                error("Must provide arguments of type FreeVariable to constructor")
            end
            fv_array[i] = free_variables[i]
        end
        new{Base.length(free_variables)}(fv_array)
    end
end

""" 
    length(::XVector)

Returns the number of free variables in the XVector
"""
numels(::XVector{D}) where {D} = D

"""
    getindex(::XVector, i::Int)

Returns the FreeVariable at index i in the given XVector
"""
Base.getindex(xv::XVector, i::Int) = xv.FVs[i]
function Base.getindex(xv::XVector, r::UnitRange{Int})
    outvec = Vector{FreeVariable}()
    for i in r
        push!(outvec, xv[i])
    end
    return XVector(outvec...)
end

"""
    iterate(::XVector)

Method for iterating through XVectors
"""
Base.iterate(xv::XVector, state=1) = state>numels(xv) ? nothing : (xv[state], state+1)

"""
    numels(xv::XVector)

Returns the length of the XVector if it is output as a Vector
"""
function Base.length(xv::XVector)
    els = 0
    for fv in xv
        els = els + length(fv)
    end
    return els
end

"""
    tovector(xv::XVector)

Returns the XVector as a Vector comprising all the elements of its component FreeVariables
"""
function tovector(xv::XVector)
    outvec = Vector{Float64}()
    for fv in xv
        append!(outvec,value(fv))
    end
    return outvec
end
