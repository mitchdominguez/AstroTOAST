# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                    PERIODIC ORBIT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
using LinearAlgebra
using MAT

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
            throw(AssertionError("traj is not periodic!"))
        end

        # Ensure that all segments of traj are continuous
        if !iscontinuous(traj, tol)
            throw(AssertionError("traj is not continuous!"))
        end

        # Ensure that thT_offset lies within [0, 2π]
        if !__within(thT_offset, 0, 2π)
            throw(BoundsError("thT must be within [0, 2pi]"))
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
    dimension(po::PeriodicOrbit{D}) where D

Return dimension of the dynamical model of traj
"""
dimension(po::PeriodicOrbit{D}) where {D} = D

"""
    x0(po::PeriodicOrbit)

Return the initial condition of the trajectory

NOTE that this is only equivalent to po(0) if `thT_offset`=0
"""
x0(po::PeriodicOrbit) = x0(get_traj(po))

"""
    dm(po::PeriodicOrbit)

Return the dynamical model of the trajectory
"""
dm(po::PeriodicOrbit) = dm(get_traj(po))

"""
    period(po::PeriodicOrbit)

Return the time of flight of the trajectory
"""
period(po::PeriodicOrbit) = tof(get_traj(po))

"""
    traj(po::PeriodicOrbit)

Return the Trajectory of the PO
"""
function traj(po::PeriodicOrbit)
    @warn "`traj` is being deprecated. Please use get_traj to obtain the Trajectory of a PeriodicOrbit"
    return copy(po.traj)
end

get_traj(po::PeriodicOrbit) = copy(po.traj)

"""
    get_traj(po::PeriodicOrbit, proptime::T; ndtime=true) where {T<:Real}

Return the Trajectory of the PO, corresponding to the propagation time
`proptime`. If `ndtime`, then treat `proptime` as a nondimensional time,
and if `ndtime` is false, then treat `proptime` as a longitudinal angle
"""
function get_traj(po::PeriodicOrbit, proptime::T; ndtime=true) where {T<:Real}
    # Convert ndtime to a nondimensional time if ndtime = false
    if ndtime==false
        println("Converting proptime from an angle to a nondimensional time")
        proptime = angle2time(po, proptime)
    end
    
    # Number of full revs to propagate for
    fullrevs = fld(proptime, period(po))
    
    # Create output Trajectory
    if fullrevs == 0
        return Trajectory(dm(po), x0(po), (0,proptime))
    else
        # TODO clean up this implementation
        out = get_traj(po)
        for i = 2:Int(fullrevs)
            append!(out, get_traj(po))
        end
        
        if maximum(tspan(out))>=proptime
            # println("1")
            # show(stdout, "text/plain", out)
            # println("---")
            # println(tof(out)%period(po))
            return out
        else
            # println("2")
            # show(stdout, "text/plain", out)
            # println("---")
            # println(tof(out)%period(po))
            append!(out, Trajectory(dm(po), x0(po), (0, proptime%period(po))))
            return out
        end
    end
end

"""
    get_traj(po::PeriodicOrbit, proptime::T1, starttime::T2; ndtime=true) where {T1<:Real, T2<:Real}

Return the Trajectory of the PO, corresponding to the propagation time
`proptime`. If `ndtime`, then treat `proptime` as a nondimensional time,
and if `ndtime` is false, then treat `proptime` as a longitudinal angle
"""
function get_traj(po::PeriodicOrbit, proptime::T1, starttime::T2; ndtime=true) where {T1<:Real, T2<:Real}

    # Convert ndtime to a nondimensional time if ndtime = false
    if ndtime==false
        println("Converting proptime and starttime from angles to nondimensional time")
        proptime = angle2time(po, proptime)
        starttime = angle2time(po, starttime)
    end

    # Convert starttime to a time within the tspan of the periodic orbit
    st = __local_time(po, starttime)
    timeleftinpo = period(po)-st
    adjustedproptime = proptime-timeleftinpo

    # Return a Trajectory, starting at starttime, but respecting the segments
    # that exist within the original periodic orbit trajectory
    potraj = get_traj(po)
    st_ind = 0
    i = 1
    seg_not_found = true

    # Find which segment within potraj that starttime is within
    while i <= length(potraj) && seg_not_found
        if __within(st, solvec(potraj)[i].t[begin], solvec(potraj)[i].t[end])
            st_ind = i
            seg_not_found = false
        end
        i+=1
    end

    if st == get_tf(solvec(potraj)[st_ind])
        st_ind+=1
    end

    # Construct Trajectory
    outtraj = Trajectory(dm(po), po(starttime; ndtime=true), [st, get_tf(solvec(potraj)[st_ind])].-st)
    # println(tof(outtraj))
    # println(timeleftinpo)
    # println(starttime)
    # println(adjustedproptime)
    # println(proptime)

    if st_ind < length(potraj)
        for i = st_ind+1:length(potraj)
            # append!(outtraj, Trajectory(dm(po), solvec(potraj)[i].u[begin], sort(tspan(solvec(potraj)[i]))))
            append!(outtraj, Trajectory(dm(po), get_x0(potraj[i]), sort(tspan(solvec(potraj)[i]))))
        end
    end

    # If adjustedproptime > 0, then get the rest of the requested trajectory
    # temp = get_traj(po, adjustedproptime; ndtime=true)
    # println(outtraj(maximum(tspan(outtraj))))
    # println(temp(0))
    # println(temp(maximum(tspan(temp))))
    # println(adjustedproptime%period(po))
    append!(outtraj, get_traj(po, adjustedproptime; ndtime=true))
    return outtraj
end


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
(po::PeriodicOrbit)(T::Real; ndtime=false) = ndtime ? get_traj(po)(__local_time(po,T)) : get_traj(po)(angle2time(po, __local_longitudinal_angle(po,T)))

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
stm(po::PeriodicOrbit, θ::Real; ndtime=false) = ndtime ? stm(get_traj(po), θ) : stm(get_traj(po), angle2time(po, θ))

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
    unit_eigs(po::PeriodicOrbit, theta_T::Real=0; ϵ::Float64=1e-4, isochronous=false, d_dim=0)

Return the unit eigenvalues and eigenvectors of the periodic orbit at the
specified longitudinal angle `theta_T`

For the unity pair, the algebraic multiplicity is 2, but the geometric
multiplicity is 1. There is only one eigenvector that gets picked up by the
eigenvalue solver. That eigenvalue, denoted `vₜ`, goes along with the dynamical
flow, as every point on a periodic orbit is itself a periodic solution.

To obtain the other eigenvector `vⱼ` which corresponds to finding the nearest
periodic solution with a different energy, the generalized eigenvector problem
must be solved `(Φ(T,0) - I)vⱼ=̇x`, where the RHS is equivalent to the time
derivative of the phase.

Setting keyword argument `isochronous=true` provides the isochronous correspondence
of the generalized eigenvector. This provides the linear prediction of a location
on a nearby PO after `theta_T` longitudinal angle of the reference orbit has passed. 
Conversely, setting `isochronous=false` provides the normal correspondence of the
generalized eigenvector at that point, essentially solving for the generalized 
eigenvector if the monodromy matrix was initially calculated at that point. Note
that this process introduces a non-negligible amount of numerical error into
the results, so be careful not to set `theta_T` too large. 

Keyword argument `d_dim` is the dimensional length (in km) of the position
component of the unity eigenvector in the direction of the flow. Setting
`d_dim` means that no length scaling occurs. Nonzero values of `d_dim` 
scale the ordinary eigenvector at `theta_T=0` to have position length `d_dim`. 
"""
function unit_eigs(po::PeriodicOrbit, theta_T::Real=0; ϵ::Float64=1e-4, isochronous=false, d_dim=0, sigma_tol=1e-10)
    # if theta_T != 0
        # @warn "unit_eigs can only return the eigenvector in the family direction if theta_T = 0!!!"
    # end
    M = monodromy(po)

    if isochronous
        # obtain isochronous correspondence of generalized eigenvector
        U,S,V = svd(M-I)
        @assert abs(S[end]) <= sigma_tol "smallest singular value cannot be approximated as zero!!!"
        S[end] = 0
        MI = U*diagm(S)*V'
        vᵥ = subspace_stepoff(normalize(dm(po)(x0(po))), d_dim, dm(po))
        vⱼ = pinv(MI)*vᵥ

        return (eigvals(po)[classify_eigs(po,ϵ)[3]], map(x->stm(po, theta_T)*x, [vᵥ, vⱼ]))
    else
        # Obtain normal correspondence of generalized eigenvector
        phi10 = stm(po, theta_T) # STM from t0 to provided theta_T

        U,S,V = svd(phi10*M*inv(phi10) - I)
        @assert abs(S[end]) <= sigma_tol "smallest singular value ($(abs(S[end]))) cannot be approximated as zero!!!"
        S[end] = 0
        MI = U*diagm(S)*V'
        vᵥ = phi10*subspace_stepoff(normalize(dm(po)(x0(po))), d_dim, dm(po))
        vⱼ = pinv(MI)*vᵥ

        return (eigvals(po)[classify_eigs(po,ϵ)[3]], [vᵥ, vⱼ])

    end
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
#

"""
    to_dict(po::PeriodicOrbit, filename::String)

Convert PeriodicOrbit `po` into a Dict
"""
function to_dict(po::PeriodicOrbit)
    model = dm(po) # Dynamical model of the periodic orbit

    outdict = Dict() # Dictionary that will be used to output results

    outdict["datatype"] = "PeriodicOrbit"
    outdict["name"] = name(po)
    outdict["family"] = family(po)


    # Add model information to the dictionary
    if typeof(model) <: Cr3bpModel
        if isnothing(model.primaries)
            outdict["P1"] = ""
            outdict["P2"] = ""
            outdict["mu"] = mass_ratio(model)

        else
            outdict["P1"] = primary_bodies(model)[1].name
            outdict["P2"] = primary_bodies(model)[2].name
            outdict["mu"] = mass_ratio(model)
        end
    end

    # Get initial conditions and times of flight for each section of the
    # trajectory
    traj = get_traj(po)
    X0 = Vector{Vector{Float64}}()
    TOF = Vector{Float64}()
    for i = 1:length(traj)
        push!(X0,get_x0(traj[i]))
        push!(TOF,tof(traj[i]))
    end

    outdict["X0"] = X0
    outdict["TOF"] = TOF
    outdict["initial_conditions"] = x0(po)
    outdict["period"] = period(po)

    return outdict
end

function to_dict(povec::Vector{PeriodicOrbit})
    outdict = Dict("data"=>[])
    for po in povec
        push!(outdict["data"], to_dict(po))
    end
    
    return outdict
end

"""
    to_mat(po::PeriodicOrbit, filename::String)

Save PeriodicOrbit `po` to a .mat file 
"""
function to_mat(po::PeriodicOrbit, filename::String)
    outdict = to_dict(po)
    endswith(filename, ".mat") ? nothing : filename = filename*".mat"
    matwrite(filename,Dict("data"=>[outdict]))
    println("Saved to $(filename)")
end

function to_mat(povec::Vector{PeriodicOrbit}, filename::String)
    outdict = to_dict(povec)
    endswith(filename, ".mat") ? nothing : filename = filename*".mat"
    matwrite(filename,outdict)
    println("Saved to $(filename)")
end

#TODO write a test case for saving and loading PeriodicOrbits
function targetcontinuity(model, _X0, _TOF, _name, _family)

    ## Create XVector
    states = Vector{FreeVariable}()
    times = Vector{FreeVariable}()
    for i = 1:length(_X0)
        if i == 1
            # Anchor states to the initial guess
            push!(states, FreeVariable("x$(i)", _X0[i], includeinds=[1,3,4,5,6]))
        else
            push!(states, FreeVariable("x$(i)", _X0[i]))
        end
        push!(times, FreeVariable("T$(i)", _TOF[i]))
    end

    xv = XVector(states..., times...)

    ## Create continuity constraints
    cc = Vector{ContinuityConstraint}()
    for i = 1:length(states)-1
        push!(cc, ContinuityConstraint(states[i], states[i+1], times[i], model))
    end
    push!(cc, ContinuityConstraint(states[end], states[1], times[end], model, includeinds=[1,2,3,4,6]))

    # Define FXVector
    fx = FXVector(cc...) # FX vector for full X, rm in FX

    # Define Targeter
    maxiter = 20
    tol = 1e-12
    targ = Targeter(xv, fx, maxiter, tol);

    # po = PeriodicOrbit(model, states, times; name = _name, family = _family)
    # po = PeriodicOrbit(model, states, times)#; name = _name, family = _family)
    # return po, targ
    try 
        # Target fx rc version
        Xhist, err = target(targ, debug=false);

        po = PeriodicOrbit(model, states, times)#; name = _name, family = _family)
        return po, targ
    catch ee
        # if isa(ee, AssertionError)
            # println("Trying again with $(length(times)+1) patch points...")
        # end 
        println("Trying again with $(length(times)+1) patch points...")

        N_new = length(times)+1 # Number of new patch points

        if N_new >= 10

            println("Times of flight:")
            show(stdout, "text/plain",_TOF)
            println()
            println("States:")
            show(stdout, "text/plain",_X0)
            println()

            ErrorException("Tried too many times to converge!") |> throw
        end

        # Re-propagate 
        sol = solve(model, SVector(_X0[1]...), (0,sum(_TOF)))

        # Get times for new patch points based on how the integrator wants to step
        inds = Int(floor(length(sol.t)/(N_new)))*[1:(N_new)...]
        inds[end] = length(sol.t)
        newtimes = [0,sol.t[inds][1:end]...]
        newtofs = diff(newtimes)

        newtimes = [0, (sum(newtofs[1:x]) for x = 1:length(newtofs))...]

        newstates = []

        # println(newtimes-sol.t[inds][1:end])
        # Sweep through each new time to get the new state associated with it
        for i = 1:length(newtimes)-1
            # For each new time, find the index of the time from the original
            # vector _TOF that it is closest to
            minind = argmin(abs.(_TOF.-newtimes[i]))
            
            # Propagate solutions from the nearest time (in either forward or 
            # reverse time) to obtain the states at the new time
            push!(newstates, solve(model, SVector(_X0[minind]...), (_TOF[minind], newtimes[i])).u[end])
        end

        # Recursively call targetcontinuity until something converges
        return targetcontinuity(model, newstates, newtofs, _name, _family)
    end

end

"""
    getpatchpoints(model::DynamicalModel, N::Int, X0::AbstractVector, TOF::AbstractVector)
    -> newstates, newtimes

Given a DynamicalModel `model` and Vectors containing patch points `X0` and
times of flight `TOF`, recalculate a new set of `N` patch points and times of flight
"""
function getpatchpoints(model::DynamicalModel, N::Int, X0::AbstractVector, TOF::AbstractVector)
    
    # Re-propagate 
    sol = solve(model, SVector(X0[1]...), (0,sum(TOF)))

    # Get times for new patch points based on how the integrator wants to step
    inds = Int(floor(length(sol.t)/(N)))*[1:(N)...]
    inds[end] = length(sol.t)
    newtimes = [0,sol.t[inds][1:end]...]
    newtofs = diff(newtimes)

    oldtimes = [0, (sum(TOF[1:x]) for x = 1:length(TOF))...]

    newstates = []

    # println(newtimes-sol.t[inds][1:end])
    # Sweep through each new time to get the new state associated with it
    for i = 1:length(newtimes)-1
    # for i = 2:length(newtimes)-1
        # For each new time, find the index of the time from the original
        # vector TOF that it is closest to
        minind = argmin(abs.(TOF.-newtimes[i]))

        # Propagate solutions from the nearest time (in either forward or 
        # reverse time) to obtain the states at the new time
        push!(newstates, solve(model, SVector(X0[minind]...), (oldtimes[minind], newtimes[i])).u[end])
    end

    return newstates, newtofs
end

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
