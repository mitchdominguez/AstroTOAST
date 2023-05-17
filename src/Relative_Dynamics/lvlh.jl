using BlockDiagonals
using OrdinaryDiffEq
using StaticArrays
#=
lvlh.jl

Dynamical model for propagating in the LVLH frame, defined with the target state relative to P2
=#
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                           LVLH MODEL                                             #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
"""
    LVLHModel

Model type for propagating relative motion in the LVLH frame in the circular restricted three body problem.
"""
struct LVLHModel{T, DQ, P} <: AutonomousDynamicalModel{12, false}
    dynamicsmodel::Cr3bpModel
    function LVLHModel(μ::T,
                       dquants::DQ=DimensionalQuantitySet(),
                       primaries::P=nothing) where {T, DQ, P}
        cr3bpmodel = Cr3bpModel(μ, dquants, primaries)
        new{T,DQ,P}(cr3bpmodel)
    end
end

"""
    dynamics_model(m::LVLHModel)
"""
dynamics_model(m::LVLHModel) = m.dynamicsmodel

"""
    model_parameters(::Cr3bpModel)

Return static vector of the model parameters (mu)
"""
model_parameters(m::LVLHModel{T}) where {T} = SVector{1, T}(mass_ratio(m))

"""
    mass_ratio(m::Cr3bpModel)

Return the mass ratio of the Cr3bp model.
"""
mass_ratio(m::LVLHModel) = m.dynamicsmodel.μ

"""
    dimensional_quantity_set(m::Cr3bpModel)

Return the dimensional quantity set associated with the Cr3bp model
"""
dimensional_quantity_set(m::LVLHModel) = m.dynamicsmodel.dquants

"""
    primary_bodies(m::Cr3bpModel)

Return primrary bodies of the Cr3bp model
"""
primary_bodies(m::LVLHModel) = m.dynamicsmodel.primaries
primary_bodies(m::LVLHModel{T, DQ, Nothing}) where {T, DQ} =
    ArgumentError("Specified Cr3bpModel does not define primaries" *
                  "(perhaps defined only with mass ratio)") |> throw

"""
    Cr3bpModel(p1::AbstractCeleestialBody, p2::AbstractCelestialBody)

Construct a CRTPB model from the two primaries.
"""
function LVLHModel(p1::AbstractCelestialBody, p2::AbstractCelestialBody)
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

    LVLHModel(
        mu,
        DimensionalQuantitySet(mass=mstar, length=lstar, time=tstar),
        (p1, p2)
    )
end

"""
    LVLHModel(dynamicsmodel::Cr3bpModel)

Create LVLHModel from Cr3bpModel
"""
function LVLHModel(dynamicsmodel::Cr3bpModel)
    LVLHModel(primary_bodies(dynamicsmodel)...)
end

# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                          PRIMARY OPS                                             #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
"""
    primary_state(m::LVLHModel, primary::T)

Return the 6 element state vector of the primary in the CRTPB rotating frame.
"""
function primary_state(m::LVLHModel, primary::T) where {T <: Cr3bpPrimary}
    primary_state(dynamics_model(m), primary)
end

"""
    distance_to_primary(m::LVLHModel, primary::T, q)

Calculate the scalar distance from a Cr3bp primary
"""
function distance_to_primary(m::LVLHModel, primary::T, q) where {T <: Cr3bpPrimary}
    distance_to_primary(dynamics_model(m), primary, q)
end

# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                        JACOBI CONSTANT                                           #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
# """
    # jacobi_constant(m::Cr3bpModel, q::T) where {T<:SArray}

# Calculate the Jacobi constant of a state in the specific model
# """
# function jacobi_constant(m::Cr3bpModel, q::T) where {T<:SArray}
    # omega = pseudopotential(m, q)
    # v2 = q[4]^2 + q[5]^2 + q[6]^2
    # 2omega - v2
# end

# """
    # jacobi_constant(m::Cr3bpModel, q::T) where {T<:Array}

# Calculate the Jacobi constant of an Array of states
# """
# function jacobi_constant(m::Cr3bpModel, q::T) where {T<:Array}
    # # Case where single state passed in
    # if eltype(q)<:Real && length(q) == dimension(m)
        # return jacobi_constant(m, SVector(Tuple(q)))
    # end

    # # Case where an array of states is passed in
    # n = length(q)
    # if eltype(q)<:SArray
        # a = Vector{eltype(eltype(q))}(undef,n)
    # else
        # println("No method existing for finding Jacobi constant of type: $(typeof(q))")
        # throw(MethodError)
    # end
    # for i = 1:n
        # a[i] = jacobi_constant(m, q[i])
    # end
    # return a
# end

# """
    # jacobi_constant(m::Cr3bpModel, q::T) where {T<:Array}

# Calculate the Jacobi constant for states output as an ODESolution
# """
# function jacobi_constant(m::Cr3bpModel, q::T) where {T<:ODESolution}
    # if length(q.u[1])==6
        # jacobi_constant(m, q.u)
    # else
        # println("No method exists for finding Jacobi constant for ODESolution"
                # *"structs with u type: $(typeof(q.u))")
        # throw(MethodError)
    # end

# end



# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                             DYNAMICS                                             #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
"""
    LVLHModel(q::AbstractArray, p::AbstractArray, t)

Evaluate the equations of motion for the circular restricted three body problem.

For the three parameter call, `p` must be an `AbstractArray` where `p[1] = mu`.
This function is provided to allow evaluation of sensitivities with respect to the mass ratio as well
as adhering to the API required by `DifferentialEquations.jl`.`

The first parameter `q` is the 12 state represented by the target state in the CR3BP frame (elements 1-6)
and the chaser relative state in the LVLH frame (elements 7-12)
"""
# function (m::Cr3bpModel)(q::AbstractArray, p::AbstractArray, t::AbstractFloat)
function (m::LVLHModel)(q::AbstractArray, p::AbstractArray, t::Real; nlprop = true, outputAmatrix=false, targetstate::T=nothing) where {T<:Union{Nothing, Trajectory, AbstractArray}}
    mu = p[1]
    # model = dynamics_model(m)

    Lmod = LVLHModel(mu)
    model = dynamics_model(Lmod)

    # Unpack target state
    if isnothing(targetstate) && length(q) == dimension(Lmod)
        # @time begin
        # xt_BCR = SVector{6}(q[1:6]) # Faster than other one for some reason
        # xt_BCR = SVector{6}([q[1:6]...]) # Same as other one for some reason
        xt_BCR = SVector{6}(copy(q[1:6])) # Same as other one for some reason

        rc_inds = 7:9
        vc_inds = 10:12
    # end # @time

    elseif isnothing(targetstate)==false && length(q) == dimension(model)
        # @time begin
        # xt_BCR = SVector{6}([targetstate(t)...])
        xt_BCR = SVector{6}(targetstate)

        rc_inds = 1:3
        vc_inds = 4:6
    # end # @time
    else
        @error "q vector should be either length 6 or 12! Length of q passed in was $(length(q))"
    end

    # show(stdout, "text/plain", xt_BCR)


    r_M = SVector{3}(xt_BCR[1:3]-[1-mu;0;0])
    Mrdot_M = SVector{3}(xt_BCR[4:6])

    # Unpack chaser state
    rho_L = SVector{3}(q[rc_inds]) # Eventually change to incorporate multiple chasers
    Lrhodot_L = SVector{3}(q[vc_inds]) # Chaser velocity [L frame]

    # @time r_M + rho_L

    # Calculate acceleration of target -->cr3bp uses coordinates centered at barycenter
    Mtargdot_M = model(xt_BCR, p, t) # q(1:6) is already wrt barycenter
    Mrddot_M = SVector{3}(Mtargdot_M[4:6]) # Target acceleration

    # Norms of states
    r = norm(r_M) # radius of target
    rdot = (1/r)*dot(r_M,Mrdot_M) # Derivative of the norm of r -- FROM FRANZINI

    # Calculate DCM: L_C_M --> vec_L = L_C_M*vec_M
    L_C_M = rot2lvlh(r_M,Mrdot_M)

    # Calculate angular momentum, derivative of angular momentum
    h_M = cross(r_M,Mrdot_M) # Angular momentum [M frame]
    h = norm(h_M) # Norm of angular momentum
    Mhdot_M = cross(r_M,Mrddot_M) # Derivative of angular momentum [M frame]
    hdot = -dot(L_C_M*Mhdot_M,[0,1,0])

    # Position of the moon wrt earth -- CIRCULAR MODEL SPECIFIC
    r_m_e_M = [1;0;0] # r_M/E in the M frame
    Mrdot_m_e_M = [0;0;0] # Mv_M/E expressed in the M frame
    r_m_e_L = L_C_M * r_m_e_M # r_M/E in the L frame

    # Target pos in the L frame
    r_L = L_C_M*r_M

    # Calculate wMI_L
    wMI_M = [0;0;1] # ang vel of M wrt I expressed in the M frame # CIRCULAR MODEL SPECIFIC
    wMI_L = L_C_M * wMI_M # ang vel of M wrt I expressed in the L frame

    # Calculate wLM_L (from Franzini eq 31)
    wLMy = -h/r^2 # j component
    wLMz = -(r/h^2)*dot(h_M,Mrddot_M) # k component
    wLM_L = [0;wLMy;wLMz] # ang vel of L wrt M expressed in the L frame

    # Calculate wLI_L
    wLI_L = wLM_L + wMI_L

    # Calculate wMIdot_L
    wMIdot_M = zeros(3) # CIRCULAR MODEL SPECIFIC
    wMIdot_L = L_C_M * wMIdot_M 

    # Calculate wMIddot_M
    wMIddot_M = [0;0;0] # CIRCULAR MODEL SPECIFIC

    # Target jerk in the M frame
    Mrdddot_M = (-2*cross(wMI_M,Mrddot_M) 
        - 3*cross(wMIdot_M,Mrdot_M)
        - cross(wMIddot_M,r_M)
        - cross(wMIdot_M,cross(wMI_M,r_M))
        - cross(wMI_M,cross(wMIdot_M,r_M))
        - cross(wMI_M,cross(wMI_M,Mrdot_M))
        - mu*ddq(r_M)*Mrdot_M
        - (1-mu)*(ddq(r_M + r_m_e_M)*(Mrdot_M + Mrdot_m_e_M)
        - ddq(r_m_e_M)*Mrdot_m_e_M))

    # Calculate wLMdot_L (from Franzini eq 32)
    wLMydot = -(1/r)*(hdot/r + 2*rdot*wLMy)
    wLMzdot = (rdot/r - 2*hdot/h)*wLMz - (r/h^2)*dot(h_M,Mrdddot_M)
    wLMdot_L = [0;wLMydot;wLMzdot]

    # Calculate wLIdot_L
    wLIdot_L = wLMdot_L + wMIdot_L - cross(wLM_L,wMI_L) # FIXED 04/20: Added cross product term


    # Calculate Lrhoddot_L
    if nlprop
        Lrhoddot_L = (-2*cross(wLI_L,Lrhodot_L)
            - cross(wLIdot_L,rho_L)
            - cross(wLI_L,cross(wLI_L,rho_L))
            + mu*(r_L/r^3 - (r_L + rho_L)/norm(r_L+rho_L)^3)
            + (1-mu)*((r_L + r_m_e_L)/norm(r_L + r_m_e_L)^3 -
                      (r_L + rho_L + r_m_e_L)/norm(r_L + rho_L + r_m_e_L)^3))
    else
        Lrhoddot_L = (-2*cross(wLI_L,Lrhodot_L)
            - cross(wLIdot_L,rho_L)
            - cross(wLI_L,cross(wLI_L,rho_L))
            - mu*ddq(r_L)*rho_L
            - (1-mu)*ddq(r_L+r_m_e_L)*rho_L)
    end


    W = crs(wLI_L)
    Wdot = crs(wLIdot_L)
    A_rrdot = -Wdot - W*W - mu*ddq(r_L) - (1-mu)*ddq(r_L+r_m_e_L)
    A = [hcat(zeros(3,3), I(3));
         hcat(A_rrdot, -2*W)]

    if nlprop==false && outputAmatrix
        rr = A*[rho_L;Lrhodot_L];
        if (norm(Lrhoddot_L-rr[4:6])>1e-12)
            @warn "Error is $(norm(Lrhoddot_L-rr[4:6]))"
        end
    end


    # Put together qdot
    if length(q) == dimension(Lmod)
        qdot = [Mrdot_M;Mrddot_M;Lrhodot_L;Lrhoddot_L]
    else
        qdot = [Lrhodot_L;Lrhoddot_L]
    end


    if outputAmatrix
        # If we want to return the A matrix
        return (qdot, A)
    else 
        # If we are only trying to propagate
        return qdot
    end

end

(m::LVLHModel)(q::AbstractArray, p::AbstractArray; kwargs...) = m(q, p, 0.0; kwargs...)
(m::LVLHModel)(q::AbstractArray; kwargs...) = m(q, model_parameters(m); kwargs...)
(m::LVLHModel)(q::AbstractArray, t::AbstractFloat; kwargs...) = m(q, model_parameters(m), t; kwargs...)

"""
    model_eoms(m::Cr3bpModel)

Return function to evaluate the Cr3bp equations of motion
"""
model_eoms(m::LVLHModel) = m


# TODO LVLH jacobian from the linearized model!!!
"""
    model_eoms_jacobian(m::LVLHModel)

Return function to evaluate the jacobian of the Cr3bp equations of motion
"""
model_eoms_jacobian(::LVLHModel) = lvlh_jacobian

"""
    lvlh_jacobian(q, p, t)

Calculate the sensitivity of the state velocity with respect to the state in the Cr3bp
"""
function lvlh_jacobian(q, p, t)
    m_lvlh = LVLHModel(p[1])
    m_underlying = dynamics_model(m_lvlh)

    A_cr3bp = model_eoms_jacobian(m_underlying)(q, p, t)
    A_lvlh = m_lvlh(q,p,t; nlprop = false, outputAmatrix = true)[2]
    
    return BlockDiagonal([A_cr3bp, A_lvlh])
end

# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                          PROPAGATION                                             #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
"""
    solve(dm::DynamicalModel, q0, tspan; [abstol=], [reltol=], [p=], [callback=])

Solve the initial value problem for the dynamical model

Note the kwarg `separate_propagation`. When true, this propagates the target and the
chaser separately, accessing the target state through a Trajectory object (using interpolation).
When false, this propagates the entire state together in one pass. There is seemingly a 
performance increase when propagating the states together as opposed to separately.
"""
function OrdinaryDiffEq.solve(dm::LVLHModel, q0, tspan;
                  separate_propagation=false,
                  abstol=DEFAULT_ABS_TOL*10,
                  reltol=DEFAULT_REL_TOL,
                  p=model_parameters(dm),
                  callback=nothing) where {D, IAD, M}

    if separate_propagation
        targetfunc = Trajectory(dynamics_model(dm), q0[1:6], tspan)
        prob = ODEProblem{false}((qq,pp,tt)->model_eoms(dm)(qq,pp,tt; targetstate=targetfunc(tt)), SVector{6}(q0[7:12]), tspan, p)
        # prob = ODEProblem{false}((qq,pp,tt)->model_eoms(dm)(qq,pp,tt; targetstate=targetfunc), SVector{6}(q0[7:12]), tspan, p)
    else
        prob = ODEProblem{false}(model_eoms(dm), SVector{12}(q0), tspan, p)
    end
    solver = DEFAULT_SOLVER
    xc = solve(prob, solver, abstol=abstol, reltol=reltol, callback=callback)

    if separate_propagation
        return (targetfunc, xc)
    else
        return xc
    end
end

# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
#                                       UTILITY FUNCTIONS                                          #
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
function crs(w::AbstractVector)
    if length(w)!=3
        @error "Length of w must be 3!"
    end
    wcrs=[ 0        -w[3]   w[2]
           w[3]     0       -w[1]
           -w[2]    w[1]    0   ]
end

function ddq(q::T) where {T<:AbstractArray}
    qnorm = norm(q);
    soln = 1/qnorm^3*(I(3)-(3*q*q')/qnorm^2);
end



"""
    Base.show

Overload the show operator to pretty print the model to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", model::LVLHModel)
    print(io, "Model: LVLH Frame Circular Restricted 3 Body Problem\n")
    print(io, "- Bodies: $(primary_bodies(model)[1].name), $(primary_bodies(model)[2].name)\n")
    print(io, "- Mass ratio: μ = $(mass_ratio(model))\n")
    print(io, "- Dimensional length: lstar = $(dimensional_length(dimensional_quantity_set(model))) km\n")
    print(io, "- Dimensional time: tstar = $(dimensional_time(dimensional_quantity_set(model))*sec2day) days")
end
