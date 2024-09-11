# """
    # AbstractCelestialBody

# Abstract base type for all celestial bodies.
# """
# abstract type AbstractCelestialBody end

# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                  ADDITIONAL INCLUDES                                   #
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
# include("naifbody.jl") 

# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                               EPHEMERIDES ACCESS METHODS                               #
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#  These methods provide an interface to the SPICE ephemerides through convenient        #
#  methods that allow for target and observer bodies to be specified as CelestialBody    #
#  types with minimal -- next to no -- overhead.                                         #
# -------------------------------------------------------------------------------------- #
"""
    ephemeris_state(; target, epoch, frame, observer, [dim=nothing])

Use SPICE toolkit to calculate the state of a target body relative to an observer.

This method leverages the SPICE.jl package which uses cspice to calculate the state of
a body relative to another in a given frame at a given epoch.
Additionally, an optional `dim` parameter is provided that, if specified, will
nondimensionalize the output state according to the dimensional quantities it defines.

## Arguments
- `target`: target body to get state of; must provide `naif_id` method
- `epoch`: observer epoch in seconds past J2000
- `frame`: reference frame of output state vector
- `observer`: observing body (base point); must provide `naif_id` method
- `dim`: optional dimensional quantities to scale output state by

See also:
[`ephemeris_position`](@ref),
[`ephemeris_velocity`](@ref)
"""
function ephemeris_state(;target, epoch, frame, observer, dim::T=nothing) where {T}
    k = spkez(target.naif_id, epoch, string(frame), "NONE", observer.naif_id)[1]
    if T != Nothing
        dl = dimensional_length(dim)
        dv = dimensional_velocity(dim)
        @inbounds k[1] /= dl
        @inbounds k[2] /= dl
        @inbounds k[3] /= dl
        @inbounds k[4] /= dv
        @inbounds k[5] /= dv
        @inbounds k[6] /= dv
    end
    vdim = @inbounds SVector{6}(
        k[1],
        k[2],
        k[3],
        k[4],
        k[5],
        k[6]
    )
end

"""
    ephemeris_position(; target, epoch, frame, observer, [dim=nothing])

Use SPICE toolkit to calculate the position of a target body relative to an observer.

This method leverages the SPICE.jl package which uses cspice to calculate the position of
a body relative to another in a given frame at a given epoch.
Additionally, an optional `dim` parameter is provided that, if specified, will
nondimensionalize the output position according to the dimensional quantities it defines.

## Arguments
- `target`: target body to get state of; must provide `naif_id` method
- `epoch`: observer epoch in seconds past J2000
- `frame`: reference frame of output position vector
- `observer`: observing body (base point); must provide `naif_id` method
- `dim`: optional dimensional quantities to scale output position by

See also:
[`ephemeris_state`](@ref),
[`ephemeris_velocity`](@ref)
"""
function ephemeris_position(;target, epoch, frame, observer, dim::T=nothing) where {T}
    k = spkezp(target.naif_id, epoch, string(frame), "NONE", observer.naif_id)[1]
    if T == Nothing
        SVector{3}(k[1], k[2], k[3])
    else
        l = dimensional_length(dim)
        SVector{3}(k[1]/ l, k[2]/ l, k[3]/l)
    end
end

"""
    ephemeris_velocity(; target, epoch, frame, observer, [dim=nothing])

Use SPICE toolkit to calculate the velocity of a target body relative to an observer.

This method leverages the SPICE.jl package which uses cspice to calculate the velocity of
a body relative to another in a given frame at a given epoch.
Additionally, an optional `dim` parameter is provided that, if specified, will
nondimensionalize the output velocity according to the dimensional quantities it defines.

## Arguments
- `target`: target body to get state of; must provide `naif_id` method
- `epoch`: observer epoch in seconds past J2000
- `frame`: reference frame of output velocity vector
- `observer`: observing body (base point); must provide `naif_id` method
- `dim`: optional dimensional quantities to scale output velocity by

See also:
[`ephemeris_state`](@ref),
[`ephemeris_position`](@ref)
"""
function ephemeris_velocity(;target, epoch, frame, observer, dim::T=nothing) where {T}
    ephemeris_state(target=target, epoch=epoch, frame=frame, observer=observer,
        dim=dim)[4:6]
end
