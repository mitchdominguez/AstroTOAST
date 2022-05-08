#=
dim.jl

Rolfe Power
=#
"""
    DimensionalQuantitySet

Set of characteristic values identifying dimensional quantites for mass, length, and time.

From these characteristic values, corresponding values are derived for:

- velocity
- acceleration
- force
- momentum

These are accessible through the `dimensional_<quantity>` methods.

See also:
[dimensional_mass](@ref),
[dimensional_time](@ref),
[dimensional_length](@ref),
[dimensional_velocity](@ref),
[dimensional_acceleration](@ref),
[dimensional_force](@ref),
[dimensional_momentum](@ref)
"""
struct DimensionalQuantitySet{T}
    mass::T
    length::T
    time::T
    velocity::T
    acceleration::T
    force::T
    momentum::T

    function DimensionalQuantitySet(mass::T, length::T, time::T) where {T}
        if mass <= 0.0
            DomainError(mass, "characteristic mass must be postive") |> throw
        elseif length <= 0.0
            DomainError(length, "characteristic length must be postive") |> throw
        elseif time <= 0.0
            DomainError(time, "characteristic time must be postive") |> throw
        end

        velocity = length / time
        acceleration = length / time / time
        force = mass * acceleration
        momentum = mass * velocity
        new{T}(
            mass,
            length,
            time,
            velocity,
            acceleration,
            force,
            momentum
        )
    end
end

"""
    dimensional_mass(::DimensionalQuantitySet)

Retrieve the characteristic mass value of the set.
"""
dimensional_mass(dq::DimensionalQuantitySet) = dq.mass
dimensional_mass(x) = dimensional_quantity_set(x) |> dimensional_mass

"""
    dimensional_length(::DimensionalQuantitySet)

Retrieve the characteristic length value of the set.
"""
dimensional_length(dq::DimensionalQuantitySet) = dq.length
dimensional_length(x) = dimensional_quantity_set(x) |> dimensional_length

"""
    dimensional_time(::DimensionalQuantitySet)

Retrieve the characteristic time value of the set.
"""
dimensional_time(dq::DimensionalQuantitySet) = dq.time
dimensional_time(x) = dimensional_quantity_set(x) |> dimensional_time

"""
    dimensional_velocity(::DimensionalQuantitySet)

Retrieve the characteristic velocity value of the set.
"""
dimensional_velocity(dq::DimensionalQuantitySet) = dq.velocity
dimensional_velocity(x) = dimensional_quantity_set(x) |> dimensional_velocity

"""
    dimensional_acceleration(::DimensionalQuantitySet)

Retrieve the characteristic acceleration value of the set.
"""
dimensional_acceleration(dq::DimensionalQuantitySet) = dq.acceleration
dimensional_acceleration(x) = dimensional_quantity_set(x) |> dimensional_acceleration

"""
    dimensional_force(::DimensionalQuantitySet)

Retrieve the characteristic force value of the set.
"""
dimensional_force(dq::DimensionalQuantitySet) = dq.force
dimensional_force(x) = dimensional_quantity_set(x) |> dimensional_force

"""
    dimensional_momentum(::DimensionalQuantitySet)

Retrieve the characteristic momentum value of the set.
"""
dimensional_momentum(dq::DimensionalQuantitySet) = dq.momentum
dimensional_momentum(x) = dimensional_quantity_set(x) |> dimensional_momentum

"""
    DimensionalQuantitySet(;[mass=1.0], [length=1.0], [time=1.0])

Construct a dimensional quantity set using the keyword arguments `mass`, `length`,
and `time`.
Any non-specified values are set to unity.
"""
function DimensionalQuantitySet(; mass=1.0, length=1.0, time=1.0)
    DimensionalQuantitySet(mass, length, time)
end

"""
    dimensionalize_state(dq::DimensionalQuantitySet, q::AbstractArray)

Assuming that q is a six element postion-velocity vector, dimensionalize
"""
function dimensionalize_state(dq, q::T) where {T <: AbstractArray}
    if length(q) < 6
        BoundsError(q, 6)
    end
    k = similar(q)
    for i = 1:3
        @inbounds k[i] = q[i] * dimensional_length(dq)
    end
    for i = 4:6
        @inbounds k[i] = q[i] * dimensional_velocity(dq)
    end
    T(k)
end

function dimensionalize_state(dq, q::T) where {T <: SArray}
    T(
        q[1] * dimensional_length(dq),
        q[2] * dimensional_length(dq),
        q[3] * dimensional_length(dq),
        q[4] * dimensional_velocity(dq),
        q[5] * dimensional_velocity(dq),
        q[6] * dimensional_velocity(dq),
    )
end

"""
    nondimensionalize_state(dq::DimensionalQuantitySet, q::AbstractArray)

Assuming that q is a six element postion-velocity vector, nondimensionalize
"""
function nondimensionalize_state(dq, q::T) where {T <: AbstractArray}
    if length(q) < 6
        BoundsError(q, 6)
    end
    k = similar(q)
    for i = 1:3
        @inbounds k[i] = q[i] / dimensional_length(dq)
    end
    for i = 4:6
        @inbounds k[i] = q[i] / dimensional_velocity(dq)
    end
    T(k)
end

function nondimensionalize_state(dq, q::T) where {T <: SArray}
    T(
        q[1] / dimensional_length(dq),
        q[2] / dimensional_length(dq),
        q[3] / dimensional_length(dq),
        q[4] / dimensional_velocity(dq),
        q[5] / dimensional_velocity(dq),
        q[6] / dimensional_velocity(dq),
    )
end

