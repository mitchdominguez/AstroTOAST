# Functions that are useful for the generation of manifolds

"""
    subpace_stepoff(v::AbstractVector, d::Real, model::DynamicalModel, theta::Real=0)

Step a distance of d dimensional position units in the subspace of the given eigenvector v
"""
function subspace_stepoff(v::AbstractVector, d_dim::Real, model::Cr3bpModel, theta::Real=0)
    # Normalize the eigenvector by the position
    v_norm = v./norm(v[1:3])
    
    # Convert d into ndim units
    d = d_dim/dimensional_length(model)

    # Stepoff vector
    u0 = d*(real(v_norm)*cos(theta) - imag(v_norm)*sin(theta))
    
    return u0
end

"""
    stable_manifold(po::PeriodicOrbit, theta_T::Real, d_dim::Real, proptime::Real)

For each stable mode of the periodic orbit, propagate the stable manifold of the periodic
orbit, starting at longitudinal angle theta_T on the PO, stepping off d_dim dimensional 
position units, and propagating for proptime nondimensional time.
"""
function stable_manifold(po::PeriodicOrbit, theta_T::Real, d_dim::Real, proptime::Real; ϵ=1e-4)
    @warn "CHECK THAT THE BACKWARDS PROPAGATION IS CORRECT"
    # Retrieve eigenvalues, eigenvectors
    lam, vee = stable_eigs(po, theta_T; ϵ) 

    # Output vector
    mans = Vector()
    q0 = 0

    # Loop through each stable mode
    for v in vee
        # Obtain stepoff
        u0 = subspace_stepoff(v, d_dim, dm(po))

        # Manifold initial state
        q0 = [po(theta_T) + u0, po(theta_T) - u0]

        # Propagate
        if typeof(dm(po)) <: Cr3bpModel
            G = diagm([1, -1, 1, -1, 1, -1])
            # tempsol = [solve(dm(po), G*q0[1], (0, proptime)), solve(dm(po), G*q0[2], (0, proptime))]
            # sol = [solve(dm(po), G*tempsol[1][end], (0, proptime)), solve(dm(po), G*tempsol[2][end], (0, proptime))]

            sol = [solve(dm(po), q0[1], (-proptime, 0)), solve(dm(po), q0[2], (-proptime, 0))]
        else
            throw(ErrorException("Negative propagation not implemented yet for this model"))
        end

        # Push to output vector
        push!(mans, sol)
    end

    return q0, mans
end

"""
    unstable_manifold(po::PeriodicOrbit, theta_T::Real, d_dim::Real, proptime::Real)

For each stable mode of the periodic orbit, propagate the stable manifold of the periodic
orbit, starting at longitudinal angle theta_T on the PO, stepping off d_dim dimensional 
position units, and propagating for proptime nondimensional time.
"""
function unstable_manifold(po::PeriodicOrbit, theta_T::Real, d_dim::Real, proptime::Real; ϵ=1e-4)
    # Retrieve eigenvalues, eigenvectors
    lam, vee = unstable_eigs(po, theta_T; ϵ) 

    # Convert time to ndtime
    tau_0 = theta_T*period(po)/(2pi)

    # Output vector
    mans = Vector()
    q0 = 0

    # Loop through each stable mode
    for v in vee
        # Obtain stepoff
        u0 = subspace_stepoff(v, d_dim, dm(po))

        # Manifold initial state
        q0 = [po(theta_T) + u0, po(theta_T) - u0]

        # Propagate
        if typeof(dm(po)) <: Cr3bpModel
            # sol = [solve(dm(po), q0[1], (0, proptime)), solve(dm(po), q0[2], (0, proptime))]
            sol = [Trajectory(dm(po), q0[1], (tau_0, tau_0+proptime)), Trajectory(dm(po), q0[2], (tau_0, tau_0+proptime))]
        else
            throw(ErrorException("Positive propagation not implemented yet for this model"))
        end

        # Push to output vector
        push!(mans, sol)
    end

    return q0, mans
end

"""
    linear_invariant_curve_2d(po::PeriodicOrbit, theta_T::Real, N::Int, d_dim::Real; ϵ=1e-4)

Generate a guess for the invariant curve of size `d_dim`, parametrized by `N` points, about a
fixed point on a periodic orbit at longitudinal angle `theta_T`
"""
function linear_invariant_curve_2d(po::PeriodicOrbit, theta_T::Real, N::Int, d_dim::Real; ϵ=1e-4)
    ths =  LinRange(0, 2pi, N)
    lam, vee = center_eigs(po, theta_T; ϵ=ϵ)

    num_modes = Int(length(lam)/2)

    u0vec = Vector{Vector{Float64}}()
    for i = 1:2:2*num_modes
        u0 = Vector{Float64}()
        for theta_rho in ths
            append!(u0, subspace_stepoff(vee[i], d_dim, dm(po), theta_rho))
        end
        push!(u0vec, u0)
    end

    return u0vec
end
