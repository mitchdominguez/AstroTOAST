"""
    solve(dm::HFEModel, q0, tspan; [abstol=], [reltol=], [p=], [callback=])

Solve the initial value problem for the high-fidelity ephemeris model
"""
function OrdinaryDiffEq.solve(dm::HFEModel, q0, tspan;
                  abstol=1e-12,
                  reltol=1e-12,
                  p=model_parameters(dm),
                  solver=DEFAULT_SOLVER,
                  callback=nothing)# where {D, IAD, M}
    prob = ODEProblem{false}(model_eoms(dm), q0, tspan, p)

    solve(prob, solver, abstol=abstol, reltol=reltol, callback=callback)
end

"""
    HFEMTangentSystem{F, JAC, N} <: Function

Structure providing equation of motion function for the augmented tangent model

    - `f` is the system EOMs
    - `J` is the Jacobian of the EOMs wrt the state, which is used to propagate the STM
    - `ws` seems to be the indices that correspond to the STM
    - `epind` is the index of the column corresponding to epoch partials
"""
struct HFEMTangentSystem{F, JAC, N} <: AbstractTangentSystem
    f::F
    j::JAC
    ws::SVector{N, Int}
    epind::Int
end


"""
    (::TangentSystem)(u, p, t)

Equations of motion for the tangent model
"""
function (tan::HFEMTangentSystem)(u, p, t)
    @inbounds s = u[:, 1]
    du = tan.f(s, p, t)
    J = tan.j(s, p, t) # A matrix 
    @inbounds dW = SMatrix{6,6,Float64}(J[:,1:(tan.ws[end]-1)] * u[:, tan.ws]) # phidot = A*phi
    @inbounds dE = J[:,1:(tan.ws[end]-1)]*u[:,tan.epind] + J[:,(tan.epind-1)]

    # # println(hcat(du, dW, dE))
    # println(typeof(du))
    # println(typeof(dW))
    # println(typeof(dE))
    # println(typeof(hcat(du, dW, dE)))

    # return nothing
    return hcat(du, dW, dE)
end

"""
        create_tangent(dm::HFEModel)

Construct a tangent model for propagating variational equations
"""
function create_tangent(dm::HFEModel)# where {D}
    eom = model_eoms(dm)
    jac = model_eoms_jacobian(dm)
    n = dimension(dm)
    ws_index = SVector{n, Int}(2:(n+1)...)
    ep_index = n+2
    tangentf = HFEMTangentSystem(eom, jac, ws_index, ep_index)
end

"""
    tangent_solve(::DynamicalModel, q0, tspan; [abstol=], [reltol=], [p=], [callback=])

Solve the initial value problem for the dynamical model propagating tangent equations
"""
function tangent_solve(dm::HFEModel, q0, tspan, Q0::M=nothing;
                       abstol=1e-12,
                       reltol=1e-12,
                       p=model_parameters(dm),
                       dense=true,
                       save_everystep=true,
                       solver = DEFAULT_SOLVER,
                       callback=nothing) where {M}
    D = dimension(dm)
    stm0 = M == Nothing ? SMatrix{D, D, eltype(q0)}(I) : Q0
    tangentf = create_tangent(dm)
    # Specify that the ode problem is not in place in type parameter
    #
    # show(stdout, "text/plain", hcat(q0, stm0, SVector{6,Float64}(zeros(6))))
    # println()
    tanprob = ODEProblem{false}(tangentf, hcat(q0, stm0, SVector{6,Float64}(zeros(6))), tspan, p)
    solve(tanprob, solver, abstol=abstol, reltol=reltol, save_everystep=save_everystep,
          dense=dense, internalnorm=_tannorm, callback=callback)

end
