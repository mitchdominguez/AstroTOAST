# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                    FREE VARIABLE
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
using StaticArrays
"""
    FreeVariable{D,T}

Type parameters include the dimensions of the model and the type of the elements
within the value vector
"""
struct FreeVariable{D,T}
    name::String
    value::Vector{T}
    # active::Bool
    removeinds::Vector{Int}
    
    # Constructor
    # function FreeVariable(name::String, value, active=true)
        # valvec = Vector{eltype(value)}()
        # append!(valvec,value)
        # new{Base.length(value), eltype(value)}(name, valvec, active) 
    # end
    function FreeVariable(name::String, value, removeinds = Vector{Int}()::Union{Int, AbstractVector{Int}})
        # Determine the dimension of the FreeVariable using value and removeinds
        rmlength = Base.length(removeinds)
        vallength = Base.length(value)

        # Check for correct bounds on removeinds
        if !all(removeinds.<=vallength) || !all(removeinds.>=1) || rmlength>vallength
            throw(BoundsError(value, removeinds))
        end
        fvdimension = vallength-rmlength # D in FreeVariable{D,T}

        # Assure that value is a vector
        valvec = Vector{eltype(value)}()
        append!(valvec,value)

        # Assure that removeinds is a vector
        rmvec = Vector{Int}()
        append!(rmvec, removeinds)

        # Construct new FreeVariable
        new{fvdimension, eltype(value)}(name, valvec, rmvec) 
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
    full_length(::FreeVariable)

Return the dimension of the free variable without the removing any elements
"""
full_length(fv::FreeVariable) = length(fv.value)

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
    removeinds(fv::FreeVariable)

Return the indices to remove from the FreeVariable
"""
removeinds(fv::FreeVariable) = copy(fv.removeinds)

"""
    active(fv::FreeVariable)

Return if the free variable is active
"""
# active(fv::FreeVariable) = fv.active
active(fv::FreeVariable{D,T}) where {D,T} = D == 0 ? false : true

"""
    value(fv::FreeVariable, rm=false)

Return the value of the free variable
"""
value(fv::FreeVariable) = deleteat!(copy(fv.value),removeinds(fv))

"""
    fullvalue(fv::FreeVariable, rm=false)

Return the value of the free variable
"""
fullvalue(fv::FreeVariable) = copy(fv.value)
# value(fv::FreeVariable, rm=false) = !rm ? copy(fv.value) : deleteat!(copy(fv.value),removeinds(fv))

"""
    tovector(fv::FreeVariable)

Return the vector value of the free variable 
"""
# tovector(fv::FreeVariable) = copy(fv.value)
tovector(fv::FreeVariable) = deleteat!(copy(fv.value),removeinds(fv))


"""
    tovector(fv::FreeVariable)

Return the vector value of the free variable without removing any elements
"""
tofullvector(fv::FreeVariable) = copy(fv.value)

"""
    tosvector(fv::FreeVariable)

Return the value of the free variable as an SVector (for increased performance with numerical integration)
"""
tosvector(fv::FreeVariable) = SVector(Tuple(tovector(fv)))

"""
    tofullsvector(fv::FreeVariable)

Return the value of the free variable as an SVector (for increased performance with numerical integration)
"""
tofullsvector(fv::FreeVariable) = SVector(Tuple(tofullvector(fv)))

"""
    getindex(::FreeVariable, i::Int)

Returns the element of the FreeVariable at index i in the given XVector,
without removing any elements.

The reason for implementing getindex in this way is so that
assignment to the entire FreeVariable can occur, and to ensure
that the only way to update a FreeVariable is to account for the 
full vector.
"""
# Base.getindex(fv::FreeVariable, i::Int) = value(fv)[i]
Base.getindex(fv::FreeVariable, i::Int) = fullvalue(fv)[i]

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
# TODO: determine if iterate needs to use the full or partial vector
# Also TODO: is iterate even necessary?
Base.iterate(fv::FreeVariable, state=1) = state>length(fv) ? nothing : (fv[state], state+1)

"""
    setindex!(fv::FreeVariable{D,T}, ind::Int, newval::T) where {D,T}

Set the value of fv at an index ind to the specified newval
"""
function Base.setindex!(fv::FreeVariable{D,T}, newval::T, ind::Int) where {D,T}
    if 1 <= ind && ind <= full_length(fv)
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
    if length(newval)!=full_length(fv)
        throw(DimensionMismatch("Update value and FreeVariable value have different dimensions"))
    end
    for i = 1:full_length(fv)
        fv[i] = newval[i]
    end
end

"""
    update(fv::FreeVariable{D,T}, newval::Vector{T}) where {D,T}

Update the values of the FreeVariable out of place
"""
function update(fv::FreeVariable{D,T}, newval::Vector{T}) where {D,T}
    if length(newval)!=full_length(fv)
        throw(DimensionMismatch("Update value and FreeVariable value have different dimensions"))
    end
    return FreeVariable(name(fv), newval)
end

"""
    Base.show

Overload the show operator to pretty print the FreeVariable to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", fv::FreeVariable)
    print(io, "Free Variable: $(name(fv))\n")
    print(io, "- Length: $(length(fv))\n")
    print(io, "- Value: $(tofullvector(fv))\n")
    print(io, "- Removed indices: $(removeinds(fv))\n")
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
        # fv_array = Vector{FreeVariable}(undef,Base.length(free_variables))
        fv_array = Vector{FreeVariable}()
        for i in eachindex(free_variables)
            if ~(typeof(free_variables[i])<:FreeVariable)
                throw(MethodError(XVector, free_variables))
            end
            # fv_array[i] = free_variables[i]
            if active(free_variables[i])
                push!(fv_array, free_variables[i])
            end
        end
        # new{Base.length(free_variables)}(fv_array)
        new{Base.length(fv_array)}(fv_array)
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
    numels(::XVector)

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
    removeinds(xv::XVector)

Output a vector containing the indices of removed elements in the full XVector
"""
function removeinds(xv::XVector)
    ind = 0
    outvec = Vector{Int}()
    for i = 1:numels(xv)
        append!(outvec, removeinds(xv[i]).+ind)
        ind += full_length(xv[i])
    end
    return outvec
end

"""
    iterate(::XVector)

Method for iterating through XVectors
"""
Base.iterate(xv::XVector, state=1) = state>numels(xv) ? nothing : (xv[state], state+1)

"""
    length(xv::XVector)

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
    full_length(xv::XVector)

Returns the length of the XVector if it is output as a Vector, without removing 
any elements
"""
function full_length(xv::XVector)
    els = 0
    for fv in xv
        els = els + full_length(fv)
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
    tofullvector(xv::XVector)

Returns the XVector as a Vector comprising all the elements of its component FreeVariables
without removing any elements
"""
function tofullvector(xv::XVector)
    outvec = Vector{Float64}()
    for fv in xv
        append!(outvec,fullvalue(fv))
    end
    return outvec
end

"""
    tosvector(xv::XVector)

Returns the XVector as an SVector comprising all the elements of its component FreeVariables
"""
function tosvector(xv::XVector)
    SVector(Tuple(tovector(xv)))
end

"""
    tofullsvector(xv::XVector)

Returns the XVector as an SVector comprising all the elements of its component FreeVariables
without removing any elements
"""
function tofullsvector(xv::XVector)
    SVector(Tuple(tofullvector(xv)))
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

Update the values of the XVector in place
"""
function update!(xv::XVector{D}, newvec::Vector{Float64}) where {D}
    if full_length(xv) == length(newvec)
        # Case that the full vector is provided in the update
        startind = 1 
        for i = 1:numels(xv)
            els = full_length(xv[i])
            update!(xv[i],newvec[startind:startind+els-1])
            startind = startind+els
        end
    elseif length(xv) == length(newvec)
        # Case that the vector with removed elements is provided in the update
        
        origvec = tofullvector(xv) # copy the original vector
        newfull = similar(origvec) # create a new vector to pass into update!
        removed = removeinds(xv) # removed inds in the XVector
        ind = 1

        for i = 1:length(newfull)
            if i in removed
                newfull[i] = origvec[i]
            else
                newfull[i] = newvec[ind]
                ind+=1
            end
        end

        update!(xv, newfull)
    else
        throw(DimensionMismatch("XVector and newvec have different number of elements"))
    end
end

"""
    update!(xv::XVector{D}, newvec::Vector{FreeVariable}) where {D}

Update the values of the FreeVariable in place
"""
function update!(xv::XVector{D}, newvec::Vector{FreeVariable}) where {D}
    update!(xv, tovector(XVector(newvec...)))
end

"""
    update!(xv::XVector{D}, xvnew::XVector{D}) where {D}

Update the values of the FreeVariable in place, given values contained within XVector xvnew
"""
function update!(xv::XVector{D}, xvnew::XVector{D}) where {D}
    update!(xv, tovector(xvnew))
end

"""
    update(xv::XVector{D}, newvec::Vector{T}) where {D,T}

Update the values of the XVector out of place
"""
function update(xv::XVector{D}, newvec::Vector{T}) where {D,T}
    xvnew = copy(xv)
    update!(xvnew, newvec)
    return xvnew
end

"""
    Base.show

Overload the show operator to pretty print the XVector to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", xv::XVector)
    print(io, "XVector:\n")
    print(io, "- Number of Free Variables: $(numels(xv))\n")
    print(io, "- Length: $(length(xv))\n")
    print(io, "---\n")
    print(io, "Free Variables:\n")
    for fv in xv
        print(io, "- Name: $(name(fv)), Length: $(length(fv))\n")
    end
end
