# Functions that are useful for the generation of manifolds
#
# TODO stable/unstable manifold propagation

"""
    subpace_stepoff(v::AbstractVector, d::Real, model::DynamicalModel, theta::Real=0)

Step a distance of d dimensional position units in the subspace of the given eigenvector v
"""
function subpace_stepoff(v::AbstractVector, d_dim::Real, model::Cr3bpModel, theta::Real=0)
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
function stable_manifold(po::PeriodicOrbit, theta_0::Real, d_dim::Real, proptime::Real)
    # Retrieve eigenvalues, eigenvectors
    lam, v = stable_eigs(po) 

    # Obtain stepoff
    # u0 = 

end
