using StaticArrays
using LinearAlgebra
#=
hfem.jl

High-Fidelity Ephemeris Model
=#

"""
    HFEModel 
"""
struct HFEModel <: NonAutonomousDynamicalModel{6, false}
    central_body::DefaultNaifBody 
    dquants::DimensionalQuantitySet
    additional_bodies::Vector{DefaultNaifBody}
    epoch::AbstractEpoch

    function HFEModel(central_body::DefaultNaifBody,
            dquants::DimensionalQuantitySet,
            additional_bodies::Vector{DefaultNaifBody},
            epoch::AbstractEpoch)

        @assert !(central_body in additional_bodies) "Central body cannot be in the list of additional bodies!"

        new(central_body, dquants, additional_bodies, epoch)
    end
end

"""
    get_central_body
"""
get_central_body(hfem::HFEModel) = hfem.central_body

"""
    dimensional_quantity_set(hfem::HFEModel)

Return the dimensional quantity set associated with the HFEM
"""
dimensional_quantity_set(hfem::HFEModel) = hfem.dquants

"""
    get_additional_bodies
"""
get_additional_bodies(hfem::HFEModel) = hfem.additional_bodies

"""
    model_parameters

Returns the universal gravitational constant [km^3/(kg.s^2)]
"""
model_parameters(::HFEModel) = @SVector [UNIVERSAL_GRAVITATIONAL_CONSTANT]

"""
    get_epoch(hfem::HFEModel)
"""
get_epoch(hfem::HFEModel) = hfem.epoch


"""
    Base.show

Overload the show operator to pretty print the model to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", model::HFEModel)
    print(io, "Model: High Fidelity Ephemeris Model\n")
    print(io, "- Central Body: $(get_central_body(model).name)\n")
    print(io, "- Additional Bodies: $([x.name for x in get_additional_bodies(model)])\n")
    # print(io, "- Mass ratio: μ = $(mass_ratio(model))\n")
    print(io, "- Dimensional length: lstar = $(dimensional_length(dimensional_quantity_set(model))) km\n")
    print(io, "- Dimensional time: tstar = $(dimensional_time(dimensional_quantity_set(model))*sec2day) days")
end
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                             DYNAMICS                                             #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
"""
    HFEModel(q::AbstractArray, p::AbstractArray, t::Real)

Evaluate the equations of motion for the high fidelity ephemeris model.

The parameter `q` is the state of the spacecraft relative to the central body, i.e, q_cs

For the three parameter call, `p` must be an `AbstractArray` where `p[1] = mu`.
This function is provided to allow evaluation of sensitivities with respect to the mass ratio as well
as adhering to the API required by `DifferentialEquations.jl`.`

In this case, `t` is the NONDIMENSIONAL time elapsed since the epoch of the HFEM

Implementation follows the equations given in the ATD mathspec
"""
function (hfem::HFEModel)(q::AbstractArray, p::AbstractArray, t::Real; debug=false)

    # @warn "Need to adjust HFEModel to take in an epoch!!!"

    G_dim = p[1] ## Universal gravitational constant -- [km^3/(kg.s^2)]

    ## Get bodies
    cb = get_central_body(hfem)
    ab = get_additional_bodies(hfem)
    if debug
        println([a.name for a in ab])
    end
    N = length(ab)

    ## Get epoch
    epoch_dim = str2et(get_epoch(hfem))

    t_dim = epoch_dim + t*dimensional_time(hfem)

    if debug
        println(et2utc(epoch_dim, :C, 14))
        println(et2utc(t_dim, :C, 14))
    end

    ## Build set of nondimensional masses
    m_c = (cb.gravitational_parameter/G_dim)/dimensional_mass(hfem)
    m_i = [bodyi.gravitational_parameter/G_dim/dimensional_mass(hfem) for bodyi in ab]

    if debug
        println(m_c)
        println(m_i)
    end

    ## Build set of state vectors from central body to additional bodies
    q_ci = [ nondimensionalize_state(hfem, ephemeris_state(;target=bodyi, epoch=t_dim,
                                                        frame="J2000",
                                                        observer=cb))
                                                        for bodyi in ab]

    ## Build set of state vectors from central
    r_ci = [x[1:3] for x in q_ci]
    r_si = [qq[1:3] - q[1:3] for qq in q_ci]
    
    r_cs = q[1:3]

    acc_dominant = -m_c/(norm(r_cs)^3)*r_cs

    acc_pert = -sum([m_i[i]*(r_ci[i]/(norm(r_ci[i])^3) - r_si[i]/(norm(r_si[i])^3)) for i=1:N])
    

    return SVector{6,Float64}(vcat(q[4:6], acc_dominant + acc_pert))

end

(m::HFEModel)(q::AbstractArray, t::AbstractFloat; debug=false) = m(q, model_parameters(m), t; debug=debug)

"""
    model_eoms(m::HFEModel)

Return function to evaluate the Cr3bp equations of motion
"""
model_eoms(m::HFEModel) = m
