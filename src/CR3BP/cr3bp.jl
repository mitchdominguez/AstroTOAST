#=
cr3bp.jl

Circular Restricted Three Body Problem (CR3BP)

Adapted from code written by Rolfe Power
=#
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                        PRIMARY BODIES                                            #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
"""
    Cr3bpPrimary

Abstract base type for Cr3bp Primaries.

Provides an "enum-like" behavior over the Cr3bp primaries.
"""
abstract type Cr3bpPrimary end

"""
    Cr3bpP1

Larger primary in the CRTPB system, at ``x = -μ``.
"""
struct Cr3bpP1 <: Cr3bpPrimary end


"""
    Cr3bpP2

Smaller primary in the CRTPB system, at ``x = 1 - μ``.
"""
struct Cr3bpP2 <: Cr3bpPrimary end
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                        MISC UTIL FUNCS                                           #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
"""
    valid_mass_ratio(μ)

Determine if the specified mass ratio is valid for a CRTPB model.

The mass ratio for a Cr3bp system is the ratio of the smaller primary's mass to the total system
mass. Consequently, it must be on the interval [0.0, 1.0].
"""
valid_mass_ratio(μ) = μ >= 0 && μ <= 1

# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                          Cr3bp MODEL                                             #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
"""
    Cr3bpModel

Model type for the circular restricted three body problem.
"""
struct Cr3bpModel{T, DQ, P} <: AutonomousDynamicalModel{6, false}
    μ::T
    dquants::DQ
    primaries::P
    function Cr3bpModel(μ::T,
                        dquants::DQ=DimensionalQuantitySet(),
                        primaries::P=nothing) where {T, DQ, P}
        if !valid_mass_ratio(μ)
            DomainError(μ, "mass ratio must be on interval [0, 1]") |> throw
        end
        new{T,DQ,P}(μ, dquants, primaries)
    end
end

"""
    model_parameters(::Cr3bpModel)

Return static vector of the model parameters (mu)
"""
model_parameters(m::Cr3bpModel{T}) where {T} = SVector{1, T}(mass_ratio(m))

"""
    mass_ratio(m::Cr3bpModel)

Return the mass ratio of the Cr3bp model.
"""
mass_ratio(m::Cr3bpModel) = m.μ

"""
    dimensional_quantity_set(m::Cr3bpModel)

Return the dimensional quantity set associated with the Cr3bp model
"""
dimensional_quantity_set(m::Cr3bpModel) = m.dquants

"""
    primary_bodies(m::Cr3bpModel)

Return primrary bodies of the Cr3bp model
"""
primary_bodies(m::Cr3bpModel) = m.primaries
primary_bodies(m::Cr3bpModel{T, DQ, Nothing}) where {T, DQ} =
    ArgumentError("Specified Cr3bpModel does not define primaries" *
                  "(perhaps defined only with mass ratio)") |> throw

"""
    Cr3bpModel(p1::AbstractCeleestialBody, p2::AbstractCelestialBody)

Construct a CRTPB model from the two primaries.
"""
function Cr3bpModel(p1::AbstractCelestialBody, p2::AbstractCelestialBody)
    if p2.parent_body != p1
        ArgumentError(
            "Incompatible primaries for Cr3bp: $(p1) and $(p2)"
        ) |> throw
    end

    gm1 = p1.gravitational_parameter
    gm2 = p2.gravitational_parameter
    mu  = gm2 / (gm1 + gm2)
    mstar = (gm1 + gm2) / UNIVERSAL_GRAVITATIONAL_CONSTANT
    lstar = p2.semimajor_axis
    tstar = sqrt(lstar^3 / (gm1 + gm2))

    Cr3bpModel(
        mu,
        DimensionalQuantitySet(mass=mstar, length=lstar, time=tstar),
        (p1, p2)
    )
end

# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                          PRIMARY OPS                                             #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
"""
    primary_state(m::Cr3bpModel, primary::Cr3bpPrimary)

Return the 6 element state vector of the primary in the CRTPB rotating frame.
"""
function primary_state(m::Cr3bpModel, primary::T) where {T <: Cr3bpPrimary}
    if T == Cr3bpP1
        SVector{6, Float64}(-mass_ratio(m), 0.0, 0.0, 0.0, 0.0, 0.0)
    elseif T == Cr3bpP2
        SVector{6, Float64}(1-mass_ratio(m), 0.0, 0.0, 0.0, 0.0, 0.0)
    else
        DomainError(primary, "primary must either be Cr3bpP1 or Cr3bpP2") |> throw
    end
end

"""
    distance_to_primary(::Cr3bpModel, primary::Cr3bpPrimary, q)

Calculate the scalar distance from a Cr3bp primary
"""
function distance_to_primary(m::Cr3bpModel, primary::T, q) where {T <: Cr3bpPrimary}
    mu = mass_ratio(m)
    if T == Cr3bpP1 
        sqrt((q[1] + mu)^2 + q[2]^2 + q[3]^2)
    elseif T == Cr3bpP2 
        sqrt((q[1] - 1 + mu)^2 + q[2]^2 + q[3]^2)
    else
        DomainError(primary, "primary must either be Cr3bpP1 or Cr3bpP2") |> throw
    end
end

# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                        PSEUDOPOTENTIAL                                           #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
"""
    pseudopotential(::Cr3bpModel, q)

Calculate the pseudopotential of a given state/position.
"""
function pseudopotential(m::Cr3bpModel, q)
    r13 = distance_to_primary(m, Cr3bpP1(), q)
    r23 = distance_to_primary(m, Cr3bpP2(), q)
    mu  = mass_ratio(m)
    (1.0 - mu) / r13 + mu / r23 + 0.5 * (q[1]^2 + q[2]^2)
end

"""
    pseudopotential_gradient(::Cr3bpModel, q)

Calculate the gradient of the pseudopotential with respect to position components.
"""
function pseudopotential_gradient(m::Cr3bpModel, q)
    r13 = distance_to_primary(m, Cr3bpP1(), q)
    r23 = distance_to_primary(m, Cr3bpP2(), q)
    r13_3 = r13^3;
    r23_3 = r23^3;
    mu  = mass_ratio(m)
    omm = 1.0 - mu
    x = q[1]
    y = q[2]
    z = q[3]
    @SVector [
       x - omm * (x + mu) / r13_3 - mu * (x - omm) / r23_3,
       y - omm * y / r13_3 - mu * y / r23_3,
       -omm * z / r13_3 - mu * z / r23_3
    ]
end

"""
    pseudopotential_jacobian(::Cr3bpModel, q)

Calculate the jacobian of the pseudopotential with respect to position components.
"""
function pseudopotential_jacobian(m::Cr3bpModel, q)
    r13 = distance_to_primary(m, Cr3bpP1(), q)
    r23 = distance_to_primary(m, Cr3bpP2(), q)
    r13_3 = r13^3;
    r23_3 = r23^3;
    r13_5 = r13_3 * r13 * r13;
    r23_5 = r23_3 * r23 * r23;
    mu  = mass_ratio(m)
    omm = 1.0 - mu
    x = q[1]
    y = q[2]
    z = q[3]

    uxx = 1 - omm/r13_3 - mu/r23_3 + 3*omm*(x+mu)^2/r13_5 + 3*mu*(x-omm)^2/r23_5;
    uyy = 1 - omm/r13_3 - mu/r23_3 + 3*omm*y^2/r13_5 + 3*mu*y^2/r23_5;
    uzz = -omm/r13_3 - mu/r23_3 + 3*omm*z^2/r13_5 + 3*mu*z^2/r23_5;

    uxy = uyx = 3*omm*(x+mu)*y/r13_5 + 3*mu*(x-omm)*y/r23_5;
    uxz = uzx = 3*omm*(x+mu)*z/r13_5 + 3*mu*(x-omm)*z/r23_5;
    uyz = uzy = 3*omm*y*z/r13_5 + 3*mu*y*z/r23_5;


    @SMatrix [
        uxx uxy uxz;
        uyx uyy uyz;
        uzx uzy uzz
    ]
end

"""
    jacobi_constant(m::Cr3bpModel, q::T) where {T<:SArray}

Calculate the Jacobi constant of a state in the specific model
"""
function jacobi_constant(m::Cr3bpModel, q::T) where {T<:SArray}
    omega = pseudopotential(m, q)
    v2 = q[4]^2 + q[5]^2 + q[6]^2
    2omega - v2
end

"""
    jacobi_constant(m::Cr3bpModel, q::T) where {T<:Array}

Calculate the Jacobi constant of an Array of states
"""
function jacobi_constant(m::Cr3bpModel, q::T) where {T<:Array}
    n = length(q)
    if eltype(q)<:SArray
        a = Vector{eltype(eltype(q))}(undef,n)
    else
        println("No method existing for finding Jacobi constant of type: $(typeof(q))")
        throw(MethodError)
    end
    for i = 1:n
        a[i] = jacobi_constant(m, q[i])
    end
    return a
end

"""
    jacobi_constant(m::Cr3bpModel, q::T) where {T<:Array}

Calculate the Jacobi constant for states output as an ODESolution
"""
function jacobi_constant(m::Cr3bpModel, q::T) where {T<:ODESolution}
    if length(q.u[1])==6
        jacobi_constant(m, q.u)
    else
        println("No method exists for finding Jacobi constant for ODESolution"
                *"structs with u type: $(typeof(q.u))")
        throw(MethodError)
    end

end



# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                             DYNAMICS                                             #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
"""
    Cr3bpModel(q::AbstractArray, p::AbstractArray, t)

Evaluate the equations of motion for the circular restricted three body problem.

For the three parameter call, `p` must be an `AbstractArray` where `p[1] = mu`.
This function is provided to allow evaluation of sensitivities with respect to the mass ratio as well
as adhering to the API required by `DifferentialEquations.jl`.`
"""
function (m::Cr3bpModel)(q::AbstractArray, p::AbstractArray, t::AbstractFloat)
    mu = p[1]
    mTemp = Cr3bpModel(mu)

    r13 = distance_to_primary(mTemp, Cr3bpP1(), q)
    r23 = distance_to_primary(mTemp, Cr3bpP2(), q)


    omm = 1.0 - mu
    r13_3 = r13^3
    r23_3 = r23^3

    @SVector [
        q[4],
        q[5],
        q[6],
        2q[5] + q[1] - omm * (q[1] + mu) / r13_3 - mu * (q[1] - omm) / r23_3,
        q[2] - 2q[4] - omm * q[2] / r13_3  - mu * q[2] / r23_3,
        -omm * q[3] / r13_3  - mu * q[3] / r23_3
    ]
end

(m::Cr3bpModel)(q::AbstractArray, p::AbstractArray) = m(q, p, 0.0)
(m::Cr3bpModel)(q::AbstractArray) = m(q, model_parameters(m))
(m::Cr3bpModel)(q::AbstractArray, t::AbstractFloat) = m(q, model_parameters(m), t)

"""
    model_eoms(m::Cr3bpModel)

Return function to evaluate the Cr3bp equations of motion
"""
model_eoms(m::Cr3bpModel) = m

"""
    model_eoms(m::Cr3bpModel)

Return function to evaluate the jacobian of the Cr3bp equations of motion
"""
model_eoms_jacobian(::Cr3bpModel) = cr3bp_jacobian

"""
    cr3bp_jacobian(q, p, t)

Calcuate the sensitivity of the state velocity with respect to the state in the Cr3bp
"""
function cr3bp_jacobian(q, p, t)
    m_local = Cr3bpModel(p[1])
    pj = pseudopotential_jacobian(m_local, q)
    @SMatrix [
            0.0     0.0     0.0  1.0 0.0 0.0;
            0.0     0.0     0.0  0.0 1.0 0.0;
            0.0     0.0     0.0  0.0 0.0 1.0;
        pj[1,1] pj[1,2] pj[1,3]  0.0 2.0 0.0;
        pj[2,1] pj[2,2] pj[2,3] -2.0 0.0 0.0;
        pj[3,1] pj[3,2] pj[3,3]  0.0 0.0 0.0
    ]
end


"""
    Base.show

Overload the show operator to pretty print the model to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", model::Cr3bpModel)
    print(io, "Model: Circular Restricted 3 Body Problem\n")
    print(io, "- Bodies: $(primary_bodies(model)[1].name), $(primary_bodies(model)[2].name)\n")
    print(io, "- Mass ratio: μ = $(mass_ratio(model))\n")
    print(io, "- Dimensional length: lstar = $(dimensional_length(dimensional_quantity_set(model))) km\n")
    print(io, "- Dimensional time: tstar = $(dimensional_time(dimensional_quantity_set(model))*sec2day) days")
end
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                      ADDITIONAL INCLUDES                                         #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
include("lagrange.jl")
