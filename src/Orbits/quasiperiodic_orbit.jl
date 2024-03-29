# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                               QUASIPERIODIC ORBIT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
using LinearAlgebra

"""
    struct QuasiPeriodicOrbit{D,N}

Object that represents quasiperiodic trajectories within a dynamical model
"""
struct QuasiPeriodicOrbit{D,N}
    ic::TrajectorySet{D,N} # Invariant curve
    T::Float64 # Stroboscopic time
    ρ::Float64 # Twist angle
    # fixedpt::Vector{Float64} # Fixed point used to generate invariant curve # TODO make this a PO or trajectory # TODO make an offset in long. angle
    fixedpt::Trajectory{D} 
    DG::Matrix{Float64} # DG matrix
    λ::Vector{ComplexF64} # Eigenvalues of DG
    V::Vector{Vector{ComplexF64}} # Eigenvectors of DG
    name::String # Name of QPO
    family::String # Family that QPO belongs to
    ##################################################
    # Define where the zero longitudinal angle of the QuasiPeriodicOrbit as defined
    # is located relative to the desired zero longitudinal angle. 
    #
    # For example, if the desired thT = 0 is located at periapsis, and the
    # QuasiPeriodicOrbit was targeted from apoapsis, then those initial
    # conditions can be passed into the constructor with thT_offset = π. 
    #
    # This prevents having to re-target orbits from where the desired thT = 0
    # actually is
    thT_offset::Float64 

    function QuasiPeriodicOrbit(ts::TrajectorySet{D,N}, rho::Float64, xstar::AbstractVector{Float64}, 
            name = "", family = "", tol=DEFAULT_CONVERGENCE_TOL, thT_offset=0) where {D,N}
        # Check that fixed point has the right dimension
        if length(xstar) != D
            throw(DimensionMismatch("dimension of xstar, and DynamicalModel are incompatible"))
        end

        # Check that N is odd
        if iseven(N)
            throw(InvalidStateException("There must be an odd number of nodes in the invariant curve", :N))
        end

        # Check that rho and T are both 1 dimensional
        if length(rho) != 1
            throw(DimensionMismatch("rho must have full length 1"))
        end
        
        ############### Check that each point on the invariant curve has the same jacobi constant ###############
        model = dm(ts)
        jc = jacobi_constant(model, x0(ts[1]))
        for traj in ts
            if abs(jacobi_constant(model, x0(traj)) - jc) > tol
                throw(ErrorException("Incompatible Jacobi constants on invariant curve"))
            end
        end
        

        ############### Check that the invariance constraint is met ###############
        # Generate FreeVariable for invariant curve
        x0vec = x0(ts)
        u0vec = similar(x0(ts))
        for i = 1:N
            u0vec[D*i-(D-1):D*i] = x0vec[D*i-(D-1):D*i] - xstar # Full state including fixed point
        end
        U0 = FreeVariable("U0", u0vec)

        # Generate FreeVariable for stroboscopic time
        T = FreeVariable("T", tof(ts))

        # Generate FreeVariable for twist angle
        rhovar = FreeVariable("rho", rho)

        # Invariance Constraint
        ic = InvarianceConstraint2D(U0, xstar, T, rhovar, dm(ts))

        if norm(evalconstraint(ic)) > tol
            println(norm(evalconstraint(ic)))
            throw(ErrorException("The given states do not correspond to a 2D torus, given the tolerance provided"))
        end

        ############### Generate TrajectorySet from -2π ≤ θ_T ≤ 2π ###############
        T = tof(ts)
        ts_backprop = TrajectorySet(dm(ts), x0(ts), T;backprop=true)

        for i = 1:N
            append!(ts_backprop[i], ts[i])
        end

        ############### Generate trajectory from fixed point ###############
        fixedpttraj = Trajectory(model, [xstar, xstar], [-tof(ts), tof(ts)])

        ############### Generate DG matrix and its eigendecomposition ###############
        DG = __dIC_du0{D*N}()(ic) + I(D*N)

        lam, vee = eigen(DG)

        V = Vector{Vector{ComplexF64}}(undef,D*N)
        for i = 1:D*N
            V[i] = vee[:,i]
        end

        # return new{D,N}(ts, tof(ts), rho, xstar, DG, lam, V, name, family, thT_offset)
        return new{D,N}(ts_backprop, T, rho, fixedpttraj, DG, lam, V, name, family, thT_offset)
    end
end

########################## OUTER CONSTRUCTORS ##########################
function QuasiPeriodicOrbit(dm::DynamicalModel, X0::Vector{Float64}, T::Float64,  
        rho::Float64, xstar::AbstractVector{Float64}; 
        name = "", family = "", tol=DEFAULT_CONVERGENCE_TOL, thT_offset=0) 
    ts = TrajectorySet(dm, X0, T)
    return QuasiPeriodicOrbit(ts, rho, xstar, name, family, tol, thT_offset)
end

########################## Basic Utilities ##########################

"""
    x0(qpo::QuasiPeriodicOrbit)

Return the initial states on the invariant curve as a Vector
"""
x0(qpo::QuasiPeriodicOrbit) = x0(qpo.ic)

"""
    xstar(qpo::QuasiPeriodicOrbit)

Return the fixed point at the initial states of the invariant curve of the QPO
"""
# xstar(qpo::QuasiPeriodicOrbit, thT=offset(qpo); ndtime=false) = ndtime ? reftraj(qpo)(thT) : reftraj(qpo)(angle2time(qpo, thT))
# xstar(qpo::QuasiPeriodicOrbit, thT=0; ndtime=false) = ndtime ? reftraj(qpo)(thT) : reftraj(qpo)(angle2time(qpo, thT))

function xstar(qpo::QuasiPeriodicOrbit, thT=0; ndtime=false)
    if ndtime
        out = reftraj(qpo)(thT)
    else
        out = reftraj(qpo)(angle2time(qpo, thT))
    end
    return out
end

function xstar(qpo::QuasiPeriodicOrbit, thT::AbstractVector; ndtime=false)
    outvec = Vector{Vector{Real}}(undef,length(thT))
    for i = 1:length(thT)
        outvec[i] = xstar(qpo, thT[i], ndtime=ndtime)
    end

    return outvec
end

"""
    reforbit(qpo::QuasiPeriodicOrbit)

Return the underlying periodic trajectory of the QPO. 

NOTE that this returns a Trajectory, and not a PeriodicOrbit
"""
reftraj(qpo::QuasiPeriodicOrbit) = qpo.fixedpt

"""
    u0(qpo::QuasiPeriodicOrbit)

Return the initial states on the invariant curve relative to the fixed point
"""
function u0(qpo::QuasiPeriodicOrbit{D,N}) where {D,N} 
    x0vec = x0(qpo)
    u0vec = similar(x0vec)
    for i = 1:N
        u0vec[D*i-(D-1):D*i] = x0vec[D*i-(D-1):D*i] - xstar(qpo) # Full state including fixed point
    end
    return u0vec
end

"""
    u2x(u::AbstractVector, xfixed::AbstractVector)

Add fixed point `xfixed` to a relative state on the qpo `u`
"""
function u2x(u::AbstractVector, xfixed::AbstractVector)
    D = length(xfixed)

    if length(u)%D != 0
        throw(ErrorException("Lengths of u and xfixed are incompatible"))
    end

    N = Int(length(u)/D)

    xvec = similar(u)
    for i = 1:N
        xvec[D*i-(D-1):D*i] = u[D*i-(D-1):D*i] + xfixed # Full state including fixed point
    end
    return xvec

end

"""
    x2u(x::AbstractVector, xfixed::AbstractVector)

Add fixed point `xfixed` to a relative state on the qpo `u`
"""
function x2u(x::AbstractVector, xfixed::AbstractVector)
    D = length(xfixed)

    if length(x)%D != 0
        throw(ErrorException("Lengths of u and xfixed are incompatible"))
    end

    N = Int(length(x)/D)

    uvec = similar(x)
    for i = 1:N
        uvec[D*i-(D-1):D*i] = x[D*i-(D-1):D*i] - xfixed # Full state including fixed point
    end
    return uvec

end

"""
    strobetime(qpo::QuasiPeriodicOrbit)

Return the stroboscopic period of the torus to be targeted
"""
strobetime(qpo::QuasiPeriodicOrbit) = qpo.T

"""
    rotationangle(qpo::QuasiPeriodicOrbit)

Return the angle by which the invariant curve will rotate after one period T
"""
rotationangle(qpo::QuasiPeriodicOrbit) = qpo.ρ

"""
    invariantcurve(qpo::QuasiPeriodicOrbit)

Return the TrajectorySet representing the invariant curve of the QPO
"""
invariantcurve(qpo::QuasiPeriodicOrbit) = qpo.ic

"""
    dm(qpo::QuasiPeriodicOrbit)

Return the dynamical model of the trajectory
"""
dm(qpo::QuasiPeriodicOrbit) = dm(invariantcurve(qpo))

"""
    dimension(qpo::QuasiPeriodicOrbit)

Return the dimension of the dynamical model of the QPO
"""
dimension(qpo::QuasiPeriodicOrbit{D,N}) where {D,N} = D

"""
    numels(qpo::QuasiPeriodicOrbit)

Return the dimension of the dynamical model of the QPO
"""
numels(qpo::QuasiPeriodicOrbit{D,N}) where {D,N} = N

"""
    offset(qpo::QuasiPeriodicOrbit)

Return the offset in longitudinal angle for `qpo` from the desired
"""
offset(qpo::QuasiPeriodicOrbit) = qpo.thT_offset

"""
    jacobi_constant(qpo::QuasiPeriodicOrbit)

Return the average jacobi constant of the states on the invariant curve
"""
function jacobi_constant(qpo::QuasiPeriodicOrbit)
    if typeof(dm(qpo)) <: Cr3bpModel 
        sum(map(x->jacobi_constant(dm(qpo), x0(x)), invariantcurve(qpo)))/length(invariantcurve(qpo)) 
    else
        throw(MethodError(jacobi_constant, qpo))
    end
end

"""
    DGmat(qpo::QuasiPeriodicOrbit)

Return the monodromy matrix of the quasiperiodic orbit
"""
DGmat(qpo::QuasiPeriodicOrbit) = qpo.DG

"""
    eigvals(qpo::QuasiPeriodicOrbit)

Return eigenvalues of the QuasiPeriodicOrbit
"""
LinearAlgebra.eigvals(qpo::QuasiPeriodicOrbit) = copy(qpo.λ)

"""
    eigvecs(qpo::QuasiPeriodicOrbit)

Return eigenvectors of the QuasiPeriodicOrbit
"""
LinearAlgebra.eigvecs(qpo::QuasiPeriodicOrbit) = copy(qpo.V)

"""
    name(qpo::QuasiPeriodicOrbit)

Print the name of the Quasiperiodic Orbit
"""
name(qpo::QuasiPeriodicOrbit) = qpo.name

"""
    family(qpo::QuasiPeriodicOrbit)

Print the family of the Quasiperiodic Orbit
"""
family(qpo::QuasiPeriodicOrbit) = qpo.family

"""
    time2angle(qpo::QuasiPeriodicOrbit)

Convert a ndim time to longitudinal angle
"""
time2angle(qpo::QuasiPeriodicOrbit, T::Real) = 2pi*T/strobetime(qpo)
time2angle(P::Real, T::Real) = 2pi*T/P

"""
    angle2time(qpo::QuasiPeriodicOrbit)

Convert a longitudinal angle to ndim time
"""
angle2time(qpo::QuasiPeriodicOrbit, th::Real) = th*strobetime(qpo)/(2pi)
angle2time(P::Real, th::Real) = th*P/(2pi)

"""
    __local_longitudinal_angle(qpo::QuasiPeriodicOrbit, th_G::Real; acceptable_range=[0,2π]::Vector)
    
Using the thT_offset, calculate what the local longitudinal angle must be
in order to output the state at the desired global longitudinal angle
"""
function __local_longitudinal_angle(qpo::QuasiPeriodicOrbit, th_G::Real; acceptable_range=[0,2π]::Vector)
    thT_offset = offset(qpo)
    if !__within(th_G, acceptable_range...)
        throw(ErrorException("th_G=$(th_G) is outside the range of acceptable values, $(acceptable_range)"))
    end

    return th_L = th_G - thT_offset
end

"""
    __local_time(qpo::QuasiPeriodicOrbit, T_G::Real; acceptable_range=[0,2π]::Vector)
    
Using the thT_offset, calculate what the local longitudinal angle must be
in order to output the state at the desired global longitudinal angle
"""
function __local_time(qpo::QuasiPeriodicOrbit, T_G::Real; acceptable_range=[0,strobetime(qpo)]::Vector)
    if !__within(T_G, acceptable_range...)
        throw(ErrorException("T_G=$(T_G) is outside the range of acceptable values, $(acceptable_range)"))
    end
    thT_offset = offset(qpo)
    T_offset = angle2time(qpo, thT_offset)
end

"""
    (qpo::QuasiPeriodicOrbit)(thT; ndtime=false)

Return the invariant curve of the quasiperiodic orbit at time T if ndtime = true, and 
returns the state at thT = 2πT/Period if ndtime = false. thT is the
longitudinal angle on the torus
"""
(qpo::QuasiPeriodicOrbit)(thT::Real; ndtime=false) = ndtime ? invariantcurve(qpo)(thT) : invariantcurve(qpo)(angle2time(qpo, thT)) ### THIS IS THE OLD IMPLEMENTATION
# (qpo::QuasiPeriodicOrbit)(thT; ndtime=false) = ndtime ? invariantcurve(qpo)(__local_time(qpo,thT)) : invariantcurve(qpo)(angle2time(qpo, __local_longitudinal_angle(qpo,thT)))


function (qpo::QuasiPeriodicOrbit)(thT::AbstractVector; ndtime=false)
    outvec = Vector{Vector{Float64}}(undef, length(thT))
    for i = 1:length(thT)
        outvec[i] = qpo(thT[i]; ndtime=ndtime)
    end
    return outvec
end

"""
    (qpo::QuasiPeriodicOrbit)(thT, thrho; ndtime=false)

For the invariant curve `ic` at longitudinal angle `thT` or time `T`, return the 
state on `ic` at latitudinal angle `thrho`
"""
function (qpo::QuasiPeriodicOrbit{dim,N})(thT::Real, thrho::Real; ndtime=false, tol = DEFAULT_CONVERGENCE_TOL) where {dim,N}
    ### NEW
    if ndtime
        # Working with ndim time
        fixedpt = xstar(qpo, thT%strobetime(qpo), ndtime=ndtime)
        uvec = x2u(qpo(thT%strobetime(qpo); ndtime=ndtime), fixedpt)
        thrho = thrho + div(thT,strobetime(qpo))*rotationangle(qpo)
    else
        # Working with angles
        fixedpt = xstar(qpo, thT%(2pi), ndtime=ndtime)
        uvec = x2u(qpo(thT%(2pi); ndtime=ndtime), fixedpt)
        thrho = thrho + div(thT,2pi)*rotationangle(qpo)
    end
    
    ### OLD
    # fixedpt = xstar(qpo, thT, ndtime=ndtime)
    # uvec = x2u(qpo(thT; ndtime=ndtime), fixedpt)


    Ut = reshape(uvec, dim, N)

    k = kvec(N)
    Dt = Matrix(transpose(Dmat(N)))


    u = Ut*Dt*Vector(vec(transpose(exp.(im*thrho*k))))

    # Remove imaginary values below a certain tolerance
    maxval, maxind = findmax(x->abs(x), imag(u))
    
    if maxval < tol
        u = real(u)
    else
        println(maxval)
        println(maxind)
        throw(InvalidStateException("Unacceptable imaginary numerical error",:u))
    end
    return u2x(u, fixedpt)

end

function (qpo::QuasiPeriodicOrbit{dim,N})(thT::Real, thrho::AbstractVector; ndtime=false, tol = DEFAULT_CONVERGENCE_TOL) where {dim,N}
    outvec = Vector{Vector{Float64}}(undef, length(thrho))
    for i = 1:length(thrho)
        outvec[i] = qpo(thT, thrho[i]; ndtime=ndtime, tol=tol)
    end
    return outvec
end

function (qpo::QuasiPeriodicOrbit{dim,N})(thT::AbstractVector, thrho::Real; ndtime=false, tol = DEFAULT_CONVERGENCE_TOL, outputmatrix=true) where {dim,N}
    if outputmatrix
        outvec = Matrix(undef, dimension(qpo), length(thT))
        for i = 1:length(thT)
            outvec[:,i] = qpo(thT[i], thrho; ndtime=ndtime, tol=tol)
        end
    else
        outvec = Vector{Vector{Float64}}(undef, length(thT))
        for i = 1:length(thT)
            outvec[i] = qpo(thT[i], thrho; ndtime=ndtime, tol=tol)
        end
    end

    return outvec
end

function (qpo::QuasiPeriodicOrbit{dim,N})(thT::AbstractVector, thrho::AbstractVector; ndtime=false, tol = DEFAULT_CONVERGENCE_TOL, outputmatrix=true) where {dim,N}
    if outputmatrix
        outvec = Vector{Matrix{Float64}}(undef, length(thrho))
        for j = 1:length(thrho)
            outvec[j] = Matrix(undef, dimension(qpo), length(thT))
            for i = 1:length(thT)
                outvec[j][:,i] = qpo(thT[i], thrho[j]; ndtime=ndtime, tol=tol)
            end
        end
    else
        # outvec = Vector{Matrix{Float64}}(undef, length(thrho))
        outvec = Vector{Vector{Vector{Float64}}}()
        for j = 1:length(thrho)
            tempvec = Vector{Vector{Float64}}(undef, length(thT))
            for i = 1:length(thT)
                tempvec[i] = qpo(thT[i], thrho[j]; ndtime=ndtime, tol=tol)
            end
            push!(outvec, tempvec)
        end

    end
    return outvec
end
"""
    Base.show

Overload the show operator to pretty print the PeriodicOrbit to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", qpo::QuasiPeriodicOrbit{D,N}) where {D,N}
    print(io, "QuasiPeriodic Orbit: $(name(qpo)) ($(family(qpo)))\n")
    print(io, "- Dimension: $(D)\n")
    print(io, "- Stroboscopic Period: $(strobetime(qpo)) ndim = $(strobetime(qpo)*dimensional_time(dm(qpo))*sec2day) days\n")
    print(io, "- Twist Angle: $(rotationangle(qpo)) rad, $(rad2deg(rotationangle(qpo))) deg\n")
    print(io, "- Longitudinal Offset: $(offset(qpo)/π)π rad\n")
end

"""
    to_dict(qpo::QuasiPeriodicOrbit)

Convert PeriodicOrbit `po` into a Dict
"""
function to_dict(qpo::QuasiPeriodicOrbit; tol=1e-10)
    model = dm(qpo) # Dynamical model of the periodic orbit

    outdict = Dict() # Dictionary that will be used to output results

    outdict["datatype"] = "QuasiPeriodicOrbit"
    outdict["name"] = name(qpo)
    outdict["family"] = family(qpo)


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
    # traj = get_traj(po)
    # X0 = Vector{Vector{Float64}}()
    # TOF = Vector{Float64}()
    # for i = 1:length(traj)
        # push!(X0,get_x0(traj[i]))
        # push!(TOF,tof(traj[i]))
    # end

    outdict["X0"] = x0(qpo)
    outdict["TOF"] = strobetime(qpo)
    outdict["rho"] = rotationangle(qpo)
    outdict["initial_conditions"] = x0(qpo)
    outdict["strobetime"] = strobetime(qpo)
    outdict["xstar"] = Vector(xstar(qpo))
    outdict["tol"] = tol ## Needed for converging QPO again
    

    return outdict
end

function to_dict(qpovec::Vector{QuasiPeriodicOrbit}; tol=1e-10)
    outdict = Dict("data"=>[])
    for qpo in qpovec
        push!(outdict["data"], to_dict(qpo; tol=tol))
    end
    
    return outdict
end

"""
    to_mat(qpo::QuasiPeriodicOrbit, filename::String)

Save QuasiPeriodicOrbit `qpo` to a .mat file 
"""
function to_mat(qpo::QuasiPeriodicOrbit, filename::String)
    outdict = to_dict(qpo)
    endswith(filename, ".mat") ? nothing : filename = filename*".mat"
    matwrite(filename,Dict("data"=>[outdict]))
    println("Saved to $(filename)")
end

function to_mat(qpovec::Vector{QuasiPeriodicOrbit}, filename::String)
    outdict = to_dict(qpovec)
    endswith(filename, ".mat") ? nothing : filename = filename*".mat"
    matwrite(filename,outdict)
    println("Saved to $(filename)")
end
