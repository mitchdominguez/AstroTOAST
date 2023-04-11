using OrdinaryDiffEq
using StaticArrays
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                            PERPENDICULAR CROSSING CONSTRAINT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
    PerpendicularCrossingConstraint

"""
struct PerpendicularCrossingConstraint{D} <: Constraint{D}
    x1::FreeVariable
    T::FreeVariable
    tmax::Real
    xz_plane_crossing
    dm::DynamicalModel
    removeinds::Vector{Int}

    # function PerpendicularCrossingConstraint(x1, x2, T, dm, rminds = Vector{Int}()::Union{Int, AbstractVector{Int}})
    function PerpendicularCrossingConstraint(dm::DynamicalModel, x1::FreeVariable, T::FreeVariable, tmax=tofullvector(T)[1]; includeinds = Vector(1:dimension(dm)))

        if full_length(x1) != dimension(dm)
            throw(DimensionMismatch("X1 and X2 must have the same full length"))
        end

        if full_length(T) != 1
            throw(DimensionMismatch("T must be of full_length 1 (currently full_length=$(full_length(T)))"))
        end

        if active(T)==false
            throw(ErrorException("T must be active for the perpendicular crossing targeter to work"))

        end

        # Find which indices to remove. Note that because of how
        # setdiff works, it is impossible to make rminds have 
        # indices that are outside the bounds of x1 or x2
        rminds = setdiff(Vector(1:dimension(dm)), includeinds)

        rmlength = Base.length(rminds)
        vallength = full_length(x1)

        # Check for correct bounds on removeinds
        if !all(rminds.<=vallength) || !all(rminds.>=1) || rmlength>vallength
            throw(BoundsError(value, removeinds))
        end

        # Calculate dimension of continuity constraint
        #   Combine removeinds(X1), removeinds(X2), and removeinds(pcc)
        # rmvec = sort(union(removeinds(x1),removeinds(x2),rminds))
        rmvec = Vector{Int}()
        append!(rmvec, rminds)

        pcdimension = vallength-length(rmvec) # D in FreeVariable{D,T}

        ## XZ-crossing event
        condition(u,t,integrator) = u[2]
        affect!(integrator) = terminate!(integrator)
        xz_plane_crossing = ContinuousCallback(condition,affect!)

        new{pcdimension}(x1, T, tmax, xz_plane_crossing, dm, rmvec)
    end
end

"""
    full_length(pcc::PerpendicularCrossingConstraint)

Return the full length of the PerpendicularCrossingConstraint, without
removing elements
"""
full_length(pcc::PerpendicularCrossingConstraint) = full_length(pcc_x1(pcc))

"""
    x1(pcc::PerpendicularCrossingConstraint)

Return x1
"""
pcc_x1(pcc::PerpendicularCrossingConstraint) = pcc.x1

"""
    tof(pcc::PerpendicularCrossingConstraint)

Return tof
"""
pcc_tof(pcc::PerpendicularCrossingConstraint) = pcc.T
pcc_T(pcc::PerpendicularCrossingConstraint) = pcc.T

"""
    pcc_tmax(pcc::PerpendicularCrossingConstraint)

Return tmax
"""
pcc_tmax(pcc::PerpendicularCrossingConstraint) = pcc.tmax

"""
    pcctspan(pcc::PerpendicularCrossingConstraint)

Return the tspan to integrate the continuity constraint for
"""
function pcctspan(pcc::PerpendicularCrossingConstraint)
    return (0, tofullvector(pcc_T(pcc))[1])
end

"""
    xz_plane_crossing(pcc::PerpendicularCrossingConstraint)

Return the ContinuousCallback function that detects when the trajectory
crosses the XZ plane and terminates the propagation at that point
"""
xz_plane_crossing(pcc::PerpendicularCrossingConstraint) = pcc.xz_plane_crossing

"""
    model(pcc::PerpendicularCrossingConstraint)

Return model
"""
dm(pcc::PerpendicularCrossingConstraint) = pcc.dm

"""
    evalconstraint(pcc::PerpendicularCrossingConstraint)

Evaluate the continuity constraint
"""
function evalconstraint(pcc::PerpendicularCrossingConstraint)
    sol = solve(dm(pcc), tofullsvector(pcc_x1(pcc)), (0, pcc_tmax(pcc)), callback = xz_plane_crossing(pcc))
    update!(pcc_T(pcc), [sol.t[end]])

    return sol.u[end]#[[4,6]]
end


# function evalconstraint(pcc::PerpendicularCrossingConstraint, xv::AbstractVector)
function evalconstraint(pcc::PerpendicularCrossingConstraint, xv)
    x1 = SVector(xv[1:6]...)
    T = xv[7]

    sol = solve(dm(pcc), x1, (0, T), callback = xz_plane_crossing(pcc))

    return sol.u[end]
end
# """
    # evalconstraint(pcc::PerpendicularCrossingConstraint, X1::FreeVariable, X2::FreeVariable)

# Evaluate the continuity constraint
# """
# function evalconstraint(pcc::PerpendicularCrossingConstraint, X1::FreeVariable, X2::FreeVariable)
    # sol = solve(dm(pcc), tofullsvector(X1), pcctspan(pcc))
    # return sol.u[end]-tofullvector(X2)
# end

"""
    partials(pcc::PerpendicularCrossingConstraint, fv::FreeVariable)

Return the matrix of partial derivatives for the partial of the constraint with
respect to the given free variable
"""
function partials(pcc::PerpendicularCrossingConstraint, fv::FreeVariable{D,T}) where {D,T}
    if fv == pcc_x1(pcc) 
        # Partial with respect to X1(T)
        return __dpcc_dx1{D}() 
    elseif  fv == pcc_tof(pcc) && active(pcc_tof(pcc))
        # Partial with respect to T
        return __dpcc_dt{D}()
    else
        # No partial
        # return __NP{D}()
        return __NP{full_length(fv)}()
    end
end

"""
    __dpcc_dx1

Partial of the continuity constraint with respect to x1(0),
which is the STM phi(0,T)
"""
struct __dpcc_dx1{D} <: Partial{D} end
function (::__dpcc_dx1{C})(pcc::PerpendicularCrossingConstraint{R}) where {R,C}
    # println(tofullsvector(pcc_x1(pcc)))
    sol = tangent_solve(dm(pcc), tofullsvector(pcc_x1(pcc)), (0, pcc_tmax(pcc)), callback=xz_plane_crossing(pcc))

    # println("TSPAN")
    # println(pcctspan(pcc))
    # println("-----")
    # show(stdout,"text/plain",sol.u[end][:,2:end])
    # println("")
    return sol.u[end][:,2:end]
end

"""
    __dpcc_dt

Partial of the continuity constraint with respect to T,
which is the derivative of the x1 state at time T
"""
struct __dpcc_dt{D} <: Partial{D} end
function (::__dpcc_dt{C})(pcc::PerpendicularCrossingConstraint{R}) where {R,C}
    sol = solve(dm(pcc), tofullsvector(pcc_x1(pcc)), (0, pcc_tmax(pcc)), callback=xz_plane_crossing(pcc))
    return dm(pcc)((sol.u[end]))
end

"""
    Base.show

Overload the show operator to pretty print the PerpendicularCrossingConstraint to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", pcc::PerpendicularCrossingConstraint)
    print(io, "Perpendicular Crossing Constraint\n")
    print(io, "- Length: $(length(pcc))\n")
    print(io, "- X1: $(name(pcc_x1(pcc)))\n")
    print(io, "- T: $(name(pcc_tof(pcc)))\n")
    print(io, "- Tmax: $(pcc_tmax(pcc))\n")
end

# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                PERIODIC ORBIT CONSTRUCTOR FROM PERPENDICULAR CROSSING
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
    PeriodicOrbit(pcc::PerpendicularCrossingConstraint, name = "", family = "", tol=DEFAULT_CONVERGENCE_TOL)

Constructor for generating a trajectory from multiple patch points
"""
function PeriodicOrbit(pcc::PerpendicularCrossingConstraint, name = "", family = "", tol=2DEFAULT_CONVERGENCE_TOL; thT_offset=0)
    # Unpack FreeVariables
    X1 = pcc_x1(pcc)
    T = pcc_T(pcc)
    model = dm(pcc)

    # Calculate positive time and negative time trajectories
    traj = Trajectory(model, tofullsvector(X1), tofullvector(T)[1])
    trajp = Trajectory(model, tofullsvector(X1), tofullvector(T)[1])
    trajn = Trajectory(model, tofullsvector(X1), (2*tofullvector(T)[1], tofullvector(T)[1]))

    # Create full orbit
    append!(traj, trajn)

    # show(stdout, "text/plain", eigen(stm(trajp)))

    S = [zeros(3,3) I(3);-I(3) zeros(3,3)]
    Om = [0 1 0;-1 0 0;0 0 0]
    V = [I(3) zeros(3,3);-Om I(3)]
    G = Diagonal([1;-1;1;-1;1;-1])
    phi_half_P = stm(trajp)
    mat1 = [zeros(3,3) -I(3);I(3) -2*Om]
    mat2 = [-2*Om I(3);-I(3) zeros(3,3)]

    M = G*mat1*phi_half_P'*mat2*G*phi_half_P # Monodromy matrix


    return PeriodicOrbit(traj, name, family, tol; thT_offset=thT_offset, M_mat = M)
end
