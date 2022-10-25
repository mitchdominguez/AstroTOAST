# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                    PERIODIC ORBIT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
using LinearAlgebra

"""
    struct PeriodicOrbit{D}

Object that represents periodic trajectories within a dynamical model
"""
struct PeriodicOrbit{D}
    traj::Trajectory{D} # Trajectory
    M::Matrix{Float64} # Monodromy matrix
    λ::Vector{ComplexF64} # Eigenvalues
    V::Vector{Vector{ComplexF64}} # Eigenvectors
    name::String # Name of orbit
    family::String # Family that orbit belongs to
    ##################################################
    # Define where the zero longitudinal angle of the PeriodicOrbit as defined
    # is located relative to the desired zero longitudinal angle. 
    #
    # For example, if the desired thT = 0 is located at periapsis, and the
    # PeriodicOrbit was targeted from apoapsis, then those initial conditions
    # can be passed into the constructor with thT_offset = π. 
    #
    # This prevents having to re-target orbits from where the desired thT = 0
    # actually is
    thT_offset::Float64 
    


    function PeriodicOrbit(traj::Trajectory{D}, name = "", family = "", tol=DEFAULT_CONVERGENCE_TOL; thT_offset=0, M_mat = nothing) where {D}
        # Ensure that traj is periodic
        if !isperiodic(traj, tol)
            throw(ErrorException("traj is not periodic!"))
        end

        # Ensure that all segments of traj are continuous
        if !iscontinuous(traj, tol)
            throw(ErrorException("traj is not continuous!"))
        end

        # Ensure that thT_offset lies within [0, 2π]
        if !__within(thT_offset, 0, 2π)
            throw(ErrorException("thT must be within [0, 2pi]"))
        end

        if isnothing(M_mat)
            # Calculate monodromy matrix
            if thT_offset == 0
                M = stm(traj)
            else
                phi_tP_t = stm(traj) # ***

                P = tof(traj)
                T_offset = thT_offset*P/(2pi)

                t0_L = wraptoperiod(0-T_offset, P)
                q0 = traj(t0_L)

                phi_t_0 = tangent_solve(dm(traj), q0, (0, T_offset)).u[end][:,2:end] # ***

                M = phi_t_0\phi_tP_t*phi_t_0
            end
        else
            # Monodromy matrix provided
            M = M_mat
        end

        # Calculate eigenvalues and eigenvectors
        λ, vec = eigen(M)

        λ = Vector{ComplexF64}(λ)
        V = Vector{Vector{ComplexF64}}(undef,D)
        for i = 1:D
            V[i] = vec[:,i]
        end

        # Make sure the eigenvalues and eigenvectors are paired and sorted properly
        sorteigs!(V, λ, M)

        # Create new PeriodicOrbit
        new{D}(copy(traj), M, λ, V, name, family, thT_offset)
    end
end

"""
    PeriodicOrbit(dm::DynamicalModel, X0, T, name = "", family = "", tol=DEFAULT_CONVERGENCE_TOL)

Constructor for generating a trajectory from multiple patch points
"""
function PeriodicOrbit(dm::DynamicalModel, X0, T, name = "", family = "", tol=DEFAULT_CONVERGENCE_TOL; thT_offset=0)
    traj = Trajectory(dm, X0, T)
    return PeriodicOrbit(traj, name, family, tol; thT_offset=thT_offset)
end

############################
# General Utilities
############################
"""
    traj(po::PeriodicOrbit)

Return the Trajectory of the PO
"""
traj(po::PeriodicOrbit) = po.traj

"""
    dimension(po::PeriodicOrbit{D}) where D

Return dimension of the dynamical model of traj
"""
dimension(po::PeriodicOrbit{D}) where {D} = D

"""
    x0(po::PeriodicOrbit)

Return the initial condition of the trajectory

NOTE that this is only equivalent to po(0) if `thT_offset`=0
"""
x0(po::PeriodicOrbit) = x0(traj(po))

"""
    dm(po::PeriodicOrbit)

Return the dynamical model of the trajectory
"""
dm(po::PeriodicOrbit) = dm(traj(po))

"""
    period(po::PeriodicOrbit)

Return the time of flight of the trajectory
"""
period(po::PeriodicOrbit) = tof(traj(po))

"""
    offset(po::PeriodicOrbit)

Return the offset in longitudinal angle for `po` from the desired
"""
offset(po::PeriodicOrbit) = po.thT_offset

"""
    jacobi_constant(po::PeriodicOrbit)

Return the jacobi constant of the periodic orbit
"""
jacobi_constant(po::PeriodicOrbit) = typeof(dm(po)) <: Cr3bpModel ? jacobi_constant(dm(po), x0(po)) : throw(MethodError(jacobi_constant, po))

"""
    (po::PeriodicOrbit)(T; ndtime=false)

Return the state on the periodic orbit at time T if ndtime = true, and 
returns the state at θ_0 = 2πT/Period if ndtime = false. θ_0 is analagous
to the longitudinal angle on a torus or the mean anomaly in a conic orbit
"""
# (po::PeriodicOrbit)(T; ndtime=false) = ndtime ? traj(po)(T) : traj(po)(angle2time(po, T))
(po::PeriodicOrbit)(T::Real; ndtime=false) = ndtime ? traj(po)(__local_time(po,T)) : traj(po)(angle2time(po, __local_longitudinal_angle(po,T)))

function (po::PeriodicOrbit)(T::AbstractVector; ndtime=false)
    outvec = Vector{Vector{Float64}}(undef, length(T))
    for i = 1:length(T)
        outvec[i] = po(T[i]; ndtime=ndtime)
    end
    return outvec
end

"""
    monodromy(po::PeriodicOrbit)

Return the monodromy matrix of the periodic orbit
"""
monodromy(po::PeriodicOrbit) = po.M

"""
    eigvals(po::PeriodicOrbit)

Return eigenvalues of the PeriodicOrbit
"""
LinearAlgebra.eigvals(po::PeriodicOrbit) = copy(po.λ)

"""
    eigvecs(po::PeriodicOrbit)

Return eigenvectors of the PeriodicOrbit
"""
LinearAlgebra.eigvecs(po::PeriodicOrbit) = copy(po.V)

"""
    name(po::PeriodicOrbit)

Print the name of the Periodic Orbit
"""
name(po::PeriodicOrbit) = po.name

"""
    family(po::PeriodicOrbit)

Print the family of the Periodic Orbit
"""
family(po::PeriodicOrbit) = po.family

"""
    time2angle(po::PeriodicOrbit)

Convert a ndim time to longitudinal angle
"""
time2angle(po::PeriodicOrbit, T::Real) = 2pi*T/period(po)

"""
    angle2time(po::PeriodicOrbit)

Convert a longitudinal angle to ndim time
"""
angle2time(po::PeriodicOrbit, th::Real) = th*period(po)/(2pi)

"""
    stm(po::PeriodicOrbit, θ::Real)

Return the state transition matrix at longitudinal angle θ. 
Calling this function with θ=2π results in the monodromy matrix
"""
stm(po::PeriodicOrbit, θ::Real; ndtime=false) = ndtime ? stm(traj(po), θ) : stm(traj(po), angle2time(po, θ))

"""
    wrapto2pi(th)

Map angles (in radians) to the range [0, 2π]. In general, 0 will map to 0.
"""
function wrapto2pi(th)
    if th>=0 && th<=2π
        return th
    end

    # if th % 2pi == 0.0 && th > 0
        # return th_wrapped = 2π
    # elseif th%2π == 0.0 && th < 0
        # return th_wrapped = 0.0
    # end

    th_wrapped = ((th % 2π) + 2π) % 2π # theta mod 2pi
end

"""
    wraptopi(th)

Map angles (in radians) to the range [-π, π].
"""
function wraptopi(th)
    return wrapto2pi(th+pi)-pi
end

"""
    wrapto360(th)

Map angles (in degrees) to the range [0, 360]
"""
function wrapto360(th)
    if th>=0 && th<=360.0
        return th
    end

    th_wrapped = ((th % 360.0) + 360.0) % 360.0 # theta mod 2pi

end

"""
    wrapto180(th)

Map angles (in degrees) to the range [-180, 180]
"""
function wrapto180(th)
    return wrapto360(th+180.0)-180.0
end

"""
    wraptoperiod(T, P)

Map times (in ndtime) to the range [0, period]. In general, 0 will map to 0.
"""
function wraptoperiod(T, P)
    if T>=0 && T<=P
        return T
    end

    # if T % 2pi == 0.0 && T > 0
        # return T_wrapped = P
    # elseif T%P == 0.0 && T < 0
        # return T_wrapped = 0.0
    # end

    T_wrapped = ((T % P) + P) % P # Teta mod 2pi
end

"""
    __local_longitudinal_angle(po::PeriodicOrbit, th_G::Real; acceptable_range=[0,2π]::Vector)
    
Using the thT_offset, calculate what the local longitudinal angle must be
in order to output the state at the desired global longitudinal angle
"""
function __local_longitudinal_angle(po::PeriodicOrbit, th_G::Real; acceptable_range=[0,2π]::Vector)
    thT_offset = offset(po)
    # if !__within(th_G, acceptable_range...)
        # throw(ErrorException("th_G=$(th_G) is outside the range of acceptable values, $(acceptable_range)"))
    # end

    return th_L = wrapto2pi(th_G - thT_offset)
end

"""
    __local_time(po::PeriodicOrbit, T_G::Real; acceptable_range=[0,2π]::Vector)
    
Using the thT_offset, calculate what the local longitudinal angle must be
in order to output the state at the desired global longitudinal angle
"""
function __local_time(po::PeriodicOrbit, T_G::Real; acceptable_range=[0,period(po)]::Vector)
    # if !__within(T_G, acceptable_range...)
        # throw(ErrorException("T_G=$(T_G) is outside the range of acceptable values, $(acceptable_range)"))
    # end
    thT_offset = offset(po)
    T_offset = angle2time(po, thT_offset)

    return t_L = wraptoperiod(T_G - T_offset, period(po))
end

######################################################
# Stability, eigendecomposition, manifolds
######################################################

"""
    paireigs(v::Vector{Vector{ComplexF64}}, M::Matrix{Float64}, λ::Vector{ComplexF64})

Pair eigenvalues with their corresponding eigenvectors.

This function will change the value of v in place
"""
function paireigs!(v::Vector{Vector{ComplexF64}}, λ::Vector{ComplexF64}, M::Matrix{Float64})
    vcopy = copy(v)
    D = length(λ)

    err = Vector{Float64}(undef, D)

    for i = 1:D
        # Calculate error for each eigenvector
        for j = 1:length(vcopy)
            err[j] = norm((M-I(D)*λ[i])*vcopy[j])
        end
        
        # Find element of vcopy that results in the minimum error
        minval, ind = findmin(err)

        # Write eigenvector corresponding to minval to index i in v
        v[i] = vcopy[ind]

        # Remove the element from vcopy and err
        deleteat!(err, ind)
        deleteat!(vcopy, ind)
    end

    return v
end

"""
    sorteigs(v::Vector{Vector{ComplexF64}}, λ::Vector{ComplexF64}, M::Matrix{Float64})

Sort pairs of eigenvalues from the pairs with the highest magnitude to pairs with the magnitude closest to 1

This function will change the value of λ and v in place
"""
function sorteigs!(v::Vector{Vector{ComplexF64}}, λ::Vector{ComplexF64}, M::Matrix{Float64})
    # Error if odd number of eigenvalues
    if length(λ)%2 != 0
        throw(DimensionMismatch("Odd number of eigenvalues"))
    end

    # Pairs of eigenvalues
    numpairs = Int(length(λ)/2)

    # Copy λ
    lamcopy = copy(λ)


    for i = 1:2:numpairs*2
        # Sort by eigenvalue magnitude
        mags = map(x->abs(x), lamcopy)

        # Find maximum 
        maxval, maxind = findmax(mags)

        # Find reciprocal of lamcopy[ind]
        recip, recipind = findmin(map(x->abs(x), lamcopy.^(-1) .- lamcopy[maxind]))

        # Write lamcopy[maxind] and lamcopy[recipind] into elements i, i+1 in λ
        λ[i] = lamcopy[maxind]
        λ[i+1] = lamcopy[recipind]

        # Remove lamcopy[maxind] and lamcopy[recipind] from lamcopy
        deleteat!(lamcopy, Tuple([x for x in sort([maxind,recipind])]))

    end

    # Re-pair eigenvalues and eigenvectors
    paireigs!(v, λ, M)

    return (λ, v)
end

"""
    classify_eigs(po::PeriodicOrbit)

Classify eigenvalues, their associated eigenvectors, into the following categories
    - unstable (magnitude > 1 + ϵ)
    - center (magnitude == 1 (±ϵ)
    - unit (value == 1 (±ϵ)
    - stable (magnitude < 1 - ϵ)

and return the indices within eigvals(po) that they correspond to.
Each eigenvalue can only belong to one of the above categories.
"""
function classify_eigs(po::PeriodicOrbit, ϵ=1e-4::Float64)
    λ = eigvals(po)

    magmin1 = map(x->abs(x), λ) .- 1

    # Separate into unstable, center, stable
    u_inds = findall(x ->x>ϵ, magmin1) # Indices of unstable eigenvalues
    c_inds = findall(x ->x<=ϵ && x>=-ϵ, magmin1) # Indices of center eigenvalues
    s_inds = findall(x ->x<-ϵ, magmin1) # Indices of stable eigenvalues

    # Find the unit eigenvalues
    val1, unit1ind = findmin(map(x->abs(x), λ.-1))
    recip, unit2ind = findmin(map(x->abs(x), λ.^(-1) .- λ[unit1ind]))
    unit_inds = [unit1ind, unit2ind]

    # Remove indices of unit eigenvalues from c_inds
    setdiff!(c_inds, unit_inds)

    return [sort(u_inds), sort(c_inds), sort(unit_inds), sort(s_inds)]
end

"""
    unstable_eigs(po::PeriodicOrbit, theta_T::Real=0; ϵ::Float64=1e-4)

Return the unstable eigenvalues and eigenvectors of the periodic orbit
at the specified longitudinal angle theta_T
"""
function unstable_eigs(po::PeriodicOrbit, theta_T::Real=0; ϵ::Float64=1e-4)
    return (eigvals(po)[classify_eigs(po,ϵ)[1]], map(x->stm(po, theta_T)*x, eigvecs(po)[classify_eigs(po,ϵ)[1]]))
end

"""
    center_eigs(po::PeriodicOrbit, ϵ=1e-4::Float64)

Return the center eigenvalues and eigenvectors of the periodic orbit (not including unit eigs)
at the specified longitudinal angle theta_T
"""
function center_eigs(po::PeriodicOrbit, theta_T::Real=0; ϵ::Float64=1e-4)
    return (eigvals(po)[classify_eigs(po,ϵ)[2]], map(x->stm(po, theta_T)*x, eigvecs(po)[classify_eigs(po,ϵ)[2]]))
end

"""
    unit_eigs(po::PeriodicOrbit, ϵ=1e-4::Float64)

Return the unit eigenvalues and eigenvectors of the periodic orbit 
at the specified longitudinal angle theta_T
"""
function unit_eigs(po::PeriodicOrbit, theta_T::Real=0; ϵ::Float64=1e-4)
    return (eigvals(po)[classify_eigs(po,ϵ)[3]], map(x->stm(po, theta_T)*x, eigvecs(po)[classify_eigs(po,ϵ)[3]]))
end

"""
    stable_eigs(po::PeriodicOrbit, ϵ=1e-4::Float64)

Return the stable eigenvalues and eigenvectors of the periodic orbit
at the specified longitudinal angle theta_T
"""
function stable_eigs(po::PeriodicOrbit, theta_T::Real=0; ϵ::Float64=1e-4)
    return (eigvals(po)[classify_eigs(po,ϵ)[4]], map(x->stm(po, theta_T)*x, eigvecs(po)[classify_eigs(po,ϵ)[4]]))
end

"""
    stability_index(po::PeriodicOrbit)

Return the stability index of the periodic orbit
"""
function stability_index(po::PeriodicOrbit{D}) where {D}
    if D%2 != 0
        throw(DimensionMismatch("Odd number of eigenvalues"))
    end
    λ = eigvals(po)
    nu = Vector{Float64}()

    for i = 1:2:D
        push!(nu, 0.5*real(λ[i] + 1/λ[i]))
    end

    return nu, max(map(x->abs(x), nu)...)
end

# TODO periapsis, apoapsis
# TODO PO continuation
# TODO export initial conditions for each segment with TOFs

"""
    Base.show

Overload the show operator to pretty print the PeriodicOrbit to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", po::PeriodicOrbit{D}) where {D}
    print(io, "Periodic Orbit: $(name(po)) ($(family(po)))\n")
    print(io, "- Dimension: $(D)\n")
    print(io, "- Period: $(period(po)) ndim = $(period(po)*dimensional_time(dm(po))*sec2day) days\n")
    print(io, "- X0: $(x0(po))\n")
    print(io, "- Longitudinal Offset: $(offset(po)/π)π rad\n")
end
