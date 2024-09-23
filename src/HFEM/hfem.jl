using StaticArrays
using LinearAlgebra
using OrdinaryDiffEq

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
# model_parameters(::HFEModel) = @SVector [UNIVERSAL_GRAVITATIONAL_CONSTANT]
model_parameters(hfem::HFEModel) = (UNIVERSAL_GRAVITATIONAL_CONSTANT,
                                    get_central_body(hfem),
                                    get_additional_bodies(hfem),
                                    get_epoch(hfem),
                                    dimensional_quantity_set(hfem))

"""
    get_epoch(hfem::HFEModel)
"""
get_epoch(hfem::HFEModel) = hfem.epoch

"""
    get_instantaneous_EM_lstar(ndtime::Float64, hfem::HFEModel)
"""
function get_instantaneous_EM_lstar(epoch_dim::Float64, hfem::HFEModel)
    return norm(ephemeris_state(;target=Moon,
                                epoch=epoch_dim+str2et(get_epoch(hfem)),
                                frame="J2000", observer=Earth)[1:3])
end

"""
    to_ephemeris_time

Convert a nondimensional time to an ephemeris time (seconds past J2000)

Setting `variable_time=true` propagates `∫(dT/dt)dt from 0 to ndtime` to obtain
the rotating-pulsating equivalent of the dimensional epoch

Setting `variable_time=false` returns `dimensional_epoch = nondimensional_epoch*tstar`
"""
function to_ephemeris_time(ndimtime::Float64, hfem::HFEModel; variable_time=true)

    if variable_time
        func(T, p, t) = sqrt(get_instantaneous_EM_lstar(T[1],hfem)^3 /
                             (UNIVERSAL_GRAVITATIONAL_CONSTANT*dimensional_mass(hfem)))

        prob = ODEProblem{false}(func, 0.0, (0.0, ndimtime))

        return str2et(get_epoch(hfem)) + solve(prob, DEFAULT_SOLVER).u[end]
    else # Constant time
        epoch_dim = str2et(get_epoch(hfem))
        return epoch_dim + ndimtime*dimensional_time(hfem)
    end

end
function to_ephemeris_time(ndimtime::AbstractVector, hfem::HFEModel; variable_time=true) 

    if variable_time
        func(T, p, t) = sqrt(get_instantaneous_EM_lstar(T[1],hfem)^3 /
                             (UNIVERSAL_GRAVITATIONAL_CONSTANT*dimensional_mass(hfem)))

        prob = ODEProblem{false}(func, 0.0, [0.0, ndimtime[end]])

        return str2et(get_epoch(hfem)) .+ solve(prob, DEFAULT_SOLVER)(ndimtime).u

    else # constant time
        return [to_ephemeris_time(t, hfem;variable_time=false) for t in ndimtime]

    end

end


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
function (hfem::HFEModel)(q::AbstractArray, p, t::Real; debug=false)

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

"""
    model_eoms(m::HFEModel)

Return function to evaluate the jacobian of the Cr3bp equations of motion
"""
model_eoms_jacobian(::HFEModel) = hfem_jacobian

"""
    hfem_jacobian(q, p, t)

Calcuate the sensitivity of the state velocity with respect to the state in the Cr3bp

Also calculate part of the partial for the senstivity of the state wrt epoch
"""
function hfem_jacobian(q::AbstractArray, p, t::Real; debug=false)
    G_dim = p[1] ## Universal gravitational constant -- [km^3/(kg.s^2)]

    ## Get bodies
    cb = p[2]
    ab = p[3]
    epoch = p[4]
    dq = p[5]

    if debug
        println([a.name for a in ab])
    end
    N = length(ab)

    ## Get epoch
    epoch_dim = str2et(epoch)

    t_dim = epoch_dim + t*dimensional_time(dq)

    if debug
        println(et2utc(epoch_dim, :C, 14))
        println(et2utc(t_dim, :C, 14))
    end

    ## Build set of nondimensional masses
    m_c = (cb.gravitational_parameter/G_dim)/dimensional_mass(dq)
    m_i = [bodyi.gravitational_parameter/G_dim/dimensional_mass(dq) for bodyi in ab]

    if debug
        println(m_c)
        println(m_i)
    end

    ## Build set of state vectors from central body to additional bodies
    q_ci = [ nondimensionalize_state(dq, ephemeris_state(;target=bodyi, epoch=t_dim,
                                                        frame="J2000",
                                                        observer=cb))
                                                        for bodyi in ab]

    ## Build set of state vectors from central
    r_ci = [x[1:3] for x in q_ci]
    v_ci = [x[4:6] for x in q_ci]
    r_si = [qq[1:3] - q[1:3] for qq in q_ci]
    
    r_cs = q[1:3]

    inds = 1:length(r_si)

    x = 1
    y = 2
    z = 3

    A11 = m_c*(3r_cs[x]^2/(norm(r_cs)^5) - 1/(norm(r_cs)^3)) + sum([m_i[i]*(3r_si[i][x]^2/norm(r_si[i])^5 - 1/norm(r_si[i])^3) for i in inds])
    A22 = m_c*(3r_cs[y]^2/(norm(r_cs)^5) - 1/(norm(r_cs)^3)) + sum([m_i[i]*(3r_si[i][y]^2/norm(r_si[i])^5 - 1/norm(r_si[i])^3) for i in inds])
    A33 = m_c*(3r_cs[z]^2/(norm(r_cs)^5) - 1/(norm(r_cs)^3)) + sum([m_i[i]*(3r_si[i][z]^2/norm(r_si[i])^5 - 1/norm(r_si[i])^3) for i in inds])

    A12 = A21 = m_c*(3r_cs[x]*r_cs[y]/(norm(r_cs)^5)) + sum([m_i[i]*(3r_si[i][x]*r_si[i][y]/norm(r_si[i])^5) for i in inds])
    A13 = A31 = m_c*(3r_cs[x]*r_cs[z]/(norm(r_cs)^5)) + sum([m_i[i]*(3r_si[i][x]*r_si[i][z]/norm(r_si[i])^5) for i in inds])
    A23 = A32 = m_c*(3r_cs[y]*r_cs[z]/(norm(r_cs)^5)) + sum([m_i[i]*(3r_si[i][y]*r_si[i][z]/norm(r_si[i])^5) for i in inds])

    
    E = zeros(6)
    for i in inds
        B11 = m_i[i]*(1/(norm(r_si[i])^3) - 3r_si[i][x]^2/(norm(r_si[i])^5) - 1/(norm(r_ci[i])^3) + 3r_ci[i][x]^2/(norm(r_ci[i])^5))
        B22 = m_i[i]*(1/(norm(r_si[i])^3) - 3r_si[i][y]^2/(norm(r_si[i])^5) - 1/(norm(r_ci[i])^3) + 3r_ci[i][y]^2/(norm(r_ci[i])^5))
        B33 = m_i[i]*(1/(norm(r_si[i])^3) - 3r_si[i][z]^2/(norm(r_si[i])^5) - 1/(norm(r_ci[i])^3) + 3r_ci[i][z]^2/(norm(r_ci[i])^5))

        B12 = B21 = 3m_i[i]*( r_ci[i][x]*r_ci[i][y]/(norm(r_ci[i])^5) - r_si[i][x]*r_si[i][y]/(norm(r_si[i])^5) ) 
        B13 = B31 = 3m_i[i]*( r_ci[i][x]*r_ci[i][z]/(norm(r_ci[i])^5) - r_si[i][x]*r_si[i][z]/(norm(r_si[i])^5) )
        B23 = B32 = 3m_i[i]*( r_ci[i][y]*r_ci[i][z]/(norm(r_ci[i])^5) - r_si[i][y]*r_si[i][z]/(norm(r_si[i])^5) )

        B = @SMatrix [ 
                      0.0     0.0     0.0
                      0.0     0.0     0.0
                      0.0     0.0     0.0
                      B11     B12     B13
                      B21     B22     B23
                      B31     B32     B33
                     ]

        if debug
            println(ab[i].name)
            show(stdout, "text/plain", B)
            println()

            println("---")
            println(1/(norm(r_si[i])^3))
            println(3r_si[i][x]^2/(norm(r_si[i])^5))
            println(1/(norm(r_ci[i])^3))
            println(3r_ci[i][x]^2/(norm(r_ci[i])^5))
            println(r_ci[i])
        end

        E += B*v_ci[i]

    end


    A = @SMatrix [
            0.0     0.0     0.0     1.0     0.0     0.0     E[1];
            0.0     0.0     0.0     0.0     1.0     0.0     E[2];
            0.0     0.0     0.0     0.0     0.0     1.0     E[3];
            A11     A12     A13     0.0     0.0     0.0     E[4];
            A21     A22     A23     0.0     0.0     0.0     E[5];
            A31     A32     A33     0.0     0.0     0.0     E[6]
    ]

    # @warn "Include epoch partials in here as well! Make a new tangent_solve method for HFEM"
    
    return A
end


include("hfem_propagation.jl")
