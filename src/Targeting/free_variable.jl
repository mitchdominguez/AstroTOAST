# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                    FREE VARIABLE
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
"""
    FreeVariable{D,T}

Type parameters include the dimensions of the model and the type of the elements
within the value vector
"""
struct FreeVariable{D,T}
    name::String
    value::Vector{T}
    
    # Constructor
    function FreeVariable(name::String, value)
        valvec = Vector{eltype(value)}()
        append!(valvec,value)
        new{Base.length(value), eltype(value)}(name, valvec) 
    end
end

"""
    copy(fv::FreeVariable)

Create a shallow copy of a FreeVariable object
"""
function Base.copy(fv::FreeVariable)
    return FreeVariable(name(fv), value(fv))
end

"""
    length(::FreeVariable)

Return the dimension of the free variable, i.e. the number of state vector elements
"""
Base.length(::FreeVariable{D,T}) where {D,T} = D

"""
    length(::FreeVariable)

Return the dimension of the free variable, i.e. the number of state vector elements
"""
Base.eltype(::FreeVariable{D,T}) where {D,T} = T

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

function Base.getindex(fv::FreeVariable, r::UnitRange{Int})
    outvec = Vector{eltype(fv)}()
    for i in r
        push!(outvec, fv[i])
    end
    return outvec
end

"""
    iterate(::FreeVariable)

Method for iterating through FreeVariables
"""
Base.iterate(fv::FreeVariable, state=1) = state>length(fv) ? nothing : (fv[state], state+1)

"""
    setindex!(fv::FreeVariable{D,T}, ind::Int, newval::T) where {D,T}

Set the value of fv at an index ind to the specified newval
"""
function Base.setindex!(fv::FreeVariable{D,T}, newval::T, ind::Int) where {D,T}
    if 1 <= ind && ind <= length(fv)
        fv.value[ind] = newval
    else
        throw(BoundsError(fv,ind)) 
    end
end

"""
    update!(fv::FreeVariable{D,T}, newval::Vector{T}) where {D,T}

Update the values of the FreeVariable in place
"""
function update!(fv::FreeVariable{D,T}, newval::Vector{T}) where {D,T}
    if length(newval)!=length(value(fv))
        throw(DimensionMismatch("Update value and FreeVariable value have different dimensions"))
    end
    for i = 1:length(value(fv))
        fv[i] = newval[i]
    end
end

"""
    update(fv::FreeVariable{D,T}, newval::Vector{T}) where {D,T}

Update the values of the FreeVariable out of place
"""
function update(fv::FreeVariable{D,T}, newval::Vector{T}) where {D,T}
    if length(newval)!=length(value(fv))
        throw(DimensionMismatch("Update value and FreeVariable value have different dimensions"))
    end
    return FreeVariable(name(fv), newval)
end


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
                throw(MethodError(XVector, free_variables))
            end
            fv_array[i] = free_variables[i]
        end
        new{Base.length(free_variables)}(fv_array)
    end
end

"""
    copy(xv::XVector)

Create a shallow copy of a XVector object
"""
function Base.copy(xv::XVector)
    outvec = Vector{FreeVariable}()
    for i in 1:numels(xv)
        xvc = copy(xv[i])
        push!(outvec, xvc)
    end
    return XVector(outvec...)
end

""" 
    length(::XVector)

Returns the number of free variables in the XVector
"""
numels(::XVector{D}) where {D} = D

"""
    getindex(::XVector, ...)

Returns the FreeVariable at index i or in UnitRange r in the given XVector
Note that getindex returns the original FreeVariable objects in the array,
so changes to each fv elsewhere will affect the output here as well
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

"""
    Base.setindex!(xv::XVector{D}, newval::FreeVariable, ind::Int) where {D}

Set the value of fv at an index ind to the specified newval
Note that this allows changing the type of the FreeVariable
"""
function Base.setindex!(xv::XVector{D}, newval::FreeVariable, ind::Int) where {D}
    if 1 <= ind && ind <= numels(xv)
        if length(xv.FVs[ind]) == length(newval)
            xv.FVs[ind] = newval
        else
            throw(DimensionMismatch("Cannot replace a FreeVariable with a FreeVariable of different length"))
        end
    else
        throw(BoundsError(xv,ind)) 
    end
end

"""
    update!(xv::XVector{D}, newvec::Vector{Float64}) where {D}

Update the values of the FreeVariable in place
"""
function update!(xv::XVector{D}, newvec::Vector{Float64}) where {D}
    if length(xv) == length(newvec)
        startind = 1 
        for i = 1:numels(xv)
            els = length(xv[i])
            update!(xv[i],newvec[startind:startind+els-1])
            startind = startind+els
        end
    else
        throw(DimensionMismatch("XVector and newvec have different number of elements"))
    end
end

"""
    update!(xv::XVector{D}, newvec::Vector{FreeVariable}) where {D}

Update the values of the FreeVariable in place
"""
function update!(xv::XVector{D}, newvec::Vector{FreeVariable}) where {D}
    if numels(xv) == length(newvec)
        for i = 1:numels(xv)
            xv[i] = newvec[i]
        end
    else
        throw(DimensionMismatch("XVector and newvec have different number of elements"))
    end
end

"""
    update(xv::XVector{D}, newvec::Vector{T}) where {D,T}

Update the values of the FreeVariable out of place
"""
function update(xv::XVector{D}, newvec::Vector{T}) where {D,T}
    xvnew = copy(xv)
    update!(xvnew, newvec)
    return xvnew
end
