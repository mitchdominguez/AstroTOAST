using SPICE

abstract type AbstractEpoch end

struct UTCEpoch <: AbstractEpoch
    year::Int64
    month::Int64
    day::Int64
    hour::Int64
    minute::Int64
    second::Float64
end
function Base.string(utc::UTCEpoch)
    datestr = string(utc.year)*"/"*string(utc.month)*"/"*string(utc.day)*" "*string(utc.hour)*":"*string(utc.minute)*":"*string(utc.second)*" UTC"
end

struct TDBEpoch <: AbstractEpoch
    year::Int64
    month::Int64
    day::Int64
    hour::Int64
    minute::Int64
    second::Float64
end

function Base.string(tdb::TDBEpoch)
    datestr = string(tdb.year)*"/"*string(tdb.month)*"/"*string(tdb.day)*" "*string(tdb.hour)*":"*string(tdb.minute)*":"*string(tdb.second)*" TDB"
end

SPICE.str2et(epoch::AbstractEpoch) = str2et(string(epoch))


jed2et(jed::Float64) = (jed - j2000())*spd()
et2jed(et::Float64) = j2000() + et/spd()
### TODO convert between TDBEpoch and UTCEpoch 

# function TDBEpoch(utc::UTCEpoch)

# end
