# Functions that are useful for the generation of manifolds
#
# TODO stable/unstable manifold propagation

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
    u0 = d*(real(v_norm)*cosd(theta) - imag(v_norm)*sind(theta))
    
    return u0
end

"""
    stable_manifold(po::PeriodicOrbit, theta_0::Real, d_dim::Real, proptime::Real)

For each stable mode of the periodic orbit, propagate the stable manifold of the periodic
orbit, starting at longitudinal angle theta_0 on the PO, stepping off d_dim dimensional 
position units, and propagating for proptime nondimensional time.
"""
function stable_manifold(po::PeriodicOrbit, theta_0::Real, d_dim::Real, proptime::Real; ϵ=1e-4)
    # Retrieve eigenvalues, eigenvectors
    lam, vee = stable_eigs(po, theta_0; ϵ) 

    # Output vector
    mans = Vector()
    q0 = 0

    # Loop through each stable mode
    for v in vee
        # Obtain stepoff
        u0 = subspace_stepoff(v, d_dim, dm(po))

        # Manifold initial state
        q0 = po(theta_0) + u0
        println(q0)

        # Propagate
        if typeof(dm(po)) <: Cr3bpModel
            G = diagm([1, -1, 1, -1, 1, -1])
            tempsol = solve(dm(po), G*q0, (0, proptime))
            sol = solve(dm(po), G*tempsol[end], (0, proptime))
            # TODO figure out a better way to propagate in negative time
        else
            throw(ErrorException("Negative propagation not implemented yet for this model"))
        end

        # Push to output vector
        push!(mans, sol)
    end

    return q0, mans
end
