#=
model.jl

Set up dynamical model objects

Adapted from code written by Rolfe Power
=#
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                    DYNAMICAL MODEL                                     #
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
"""
    DynamicalModel{D, IAD}

Abstract base type for all dynamical models.

## Type Parameters
- `D`: dimension of the model
- `IAD`: whether or not the model jacobian is auto-differentiated
"""
abstract type DynamicalModel{D, IAD} end

"""
    dimension(::DynamicalModel)

Return the dimension of the dynamical model, i.e. the number of state vector elements
"""
dimension(::DynamicalModel{D}) where {D} = D


"""
    isautonomous(::DynamicalModel)

Return whether the model is autonomous, i.e., time-independent
"""
isautonomous(::DynamicalModel) = false

"""
    isautodiff(::DynamicalModel)

Return whether the model's jacobian must be automatically differentiated
"""
isautodiff(::DynamicalModel{D, IAD}) where {D, IAD} = IAD

"""
    model_eoms(::DynamicalModel)

Return the equations of motion for the dynamical model.

The equations of motion are passable to the `DifferentialEquations` suite.
"""
function model_eoms(dm::T) where {T <: DynamicalModel}
    MethodError(model_eoms, dm)
end

# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                              AUTONOMOUS DYANMICAL MODEL                                #
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
"""
    AutonomousDynamicalModel{D, IAD} <: DynamicalModel{D, IAD}

Abstract base type for autonomous (time independent) dynamical models

See also:
[`DynamicalModel`](@ref)
"""
abstract type AutonomousDynamicalModel{D, IAD} <: DynamicalModel{D, IAD} end

isautonomous(::AutonomousDynamicalModel) = true

# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                             JACOBIAN CREATION FOR AUTODIFF                             #
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
# This formulation is based on the excellent DynamicalSystems.jl package that can be     #
# found at the following url: https://github.com/JuliaDynamics/DynamicalSystemsBase.jl   #
# -------------------------------------------------------------------------------------- #
"""
    create_jacobian(::DynamicalModel)

Construct a jacobian function for the given dynamical model using auto-diff

Uses the automatic differentiation capabilities offered by `ForwardDiff` to return a
function that will evaluate the jacobian of the dynamical model.
"""
function create_jacobian(dm::DynamicalModel{D}) where {D}
    f = model_eoms(dm)
    if dimension(dm) == 1
        (u, p, t) -> ForwardDiff.derivative(x -> f(x, p, t), u)
    else
        (u, p, t) -> ForwardDiff.jacobian(x -> f(x, p, t), u)
    end
end

"""
    model_eoms_jacobian(::DynamicalModel{D, true})

Construct the model EOM jacobian using automatic differentiation for all IAD models.
"""
function model_eoms_jacobian(dm::DynamicalModel{D, true}) where {D}
    create_jacobian(dm)
end

"""
    solve(dm::DynamicalModel, q0, tspan; [abstol=], [reltol=], [p=], [callback=])

Solve the initial value problem for the dynamical model
"""
function OrdinaryDiffEq.solve(dm::DynamicalModel, q0, tspan;
                  abstol=DEFAULT_ABS_TOL,
                  reltol=DEFAULT_REL_TOL,
                  p=model_parameters(dm),
                  callback=nothing) where {D, IAD, M}
    prob = ODEProblem{false}(model_eoms(dm), q0, tspan, p)
    solver = DEFAULT_SOLVER
    solve(prob, solver, abstol=abstol, reltol=reltol, callback=callback)
end

# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                              TANGENT SYSTEM CONSTRUCTION                               #
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
"""
    TangentSystem{F, JAC, N} <: Function

Structure providing equation of motion function for the augmented tangent model
"""
struct TangentSystem{F, JAC, N} <: Function
    f::F
    j::JAC
    ws::SVector{N, Int}
end

"""
    (::TangentSystem)(u, p, t)

Equations of motion for the tangent model
"""
function (tan::TangentSystem)(u, p, t)
    @inbounds s = u[:, 1]
    du = tan.f(s, p, t)
    J = tan.j(s, p, t)
    @inbounds dW = J * u[:, tan.ws]
    return hcat(du, dW)
end

"""
        create_tangent(dm::DynamicalModel)

Construct a tangent model for propagating variational equations
"""
function create_tangent(dm::DynamicalModel) where {D}
    eom = model_eoms(dm)
    jac = model_eoms_jacobian(dm)
    n = dimension(dm)
    ws_index = SVector{n, Int}(2:(n+1)...)
    tangentf = TangentSystem(eom, jac, ws_index)
end

"""
    _tannorm(u::AbstractMatrix, t)

Norm function to use for tangent propagation, ignores STM
"""
function _tannorm(u::AbstractMatrix, t)
    s = size(u)[1]
    x = zero(eltype(u))
    @inbounds for i in 1:s
        x += u[i, 1]^2
    end
    return sqrt(x/length(x))
end
_tannorm(u::Real, t) = abs(u)

"""
    tangent_solve(::DynamicalModel, q0, tspan; [abstol=], [reltol=], [p=], [callback=])

Solve the initial value problem for the dynamical model propagating tangent equations
"""
function tangent_solve(dm::DynamicalModel{D, IAD}, q0, tspan, Q0::M=nothing;
                       abstol=DEFAULT_ABS_TOL,
                       reltol=DEFAULT_REL_TOL,
                       p=model_parameters(dm),
                       dense=true,
                       save_everystep=true,
                       callback=nothing) where {D, IAD, M}
    stm0 = M == Nothing ? SMatrix{D, D, eltype(q0)}(I) : Q0
    tangentf = create_tangent(dm)
    # Specify that the ode problem is not in place in type parameter
    tanprob = ODEProblem{false}(tangentf, hcat(q0, stm0), tspan, p)
    solver = DEFAULT_SOLVER
    solve(tanprob, solver, abstol=abstol, reltol=reltol, save_everystep=save_everystep,
          dense=dense, internalnorm=_tannorm, callback=callback)
end
