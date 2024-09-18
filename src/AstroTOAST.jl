module AstroTOAST

using StaticArrays
using ForwardDiff
using OrdinaryDiffEq
using DiffEqBase
using LinearAlgebra
using SPICE

include("constants.jl")

include("dim.jl")
export DimensionalQuantitySet
export dimensional_mass, dimensional_length, dimensional_time
export dimensional_velocity, dimensional_acceleration
export dimensional_force, dimensional_momentum
export dimensionalize_state, nondimensionalize_state
export dimensional_quantity_set

include("units.jl")
export m2km, km2m, au2km, km2au
export min2sec, sec2min, hr2min, min2hr, day2hr, hr2day,
       sec2hr, hr2sec, sec2day, day2sec, min2day, day2min,
       yr2day, day2yr, yr2hr, hr2yr, yr2min, min2yr, yr2sec, sec2yr

include("model.jl")
export dimension, isautonomous, isautodiff, model_parameters
export model_eoms, model_eoms_jacobian
export create_jacobian
export solve
export tangent_solve
export equilibrium_solutions

include("time.jl")
export UTCEpoch, TDBEpoch
export jed2et, et2jed

include("Body/naifbody.jl")
export Bodies

include("ephemeris.jl")
export defaultLSK, defaultSPK
export load_default_kernels, currently_loaded_kernels, clear_all_kernels, load_only_default_kernels

include("Body/body.jl")
export ephemeris_state, ephemeris_position, ephemeris_velocity

include("CR3BP/cr3bp.jl")
export Cr3bpModel
export mass_ratio, primary_bodies
export distance_to_primary, primary_state
export pseudopotential, pseudopotential_gradient, pseudopotential_jacobian
export jacobi_constant
export lyapunov_linear_ics

include("HFEM/hfem.jl")
export HFEModel
export get_central_body, get_additional_bodies, get_epoch, to_ephemeris_time

include("Targeting/free_variable.jl")
export FreeVariable
export full_length, name, fullvalue, tofullvector, value, update, active, removeinds
export XVector
export numels, tovector, tosvector, tofullsvector

include("Targeting/constraint.jl")
export FXVector

include("Targeting/targeter.jl")
export Targeter
export X, FX, DFX, getmaxiter, gettol, evalDFXMatrix, target

##################################################################
########################### ORBITS ##############################
##################################################################

include("Orbits/trajectory.jl")
export Trajectory
export x0, tspan, solvec, isperiodic, iscontinuous, stm, get_x0, get_xf, get_t0, get_tf

include("Orbits/trajectory_set.jl")
export TrajectorySet

include("Orbits/periodic_orbit.jl")
export PeriodicOrbit
export offset
export family
export period, get_traj, monodromy, lambda, paireigs!, sorteigs!
export angle2time, time2angle
export unstable_eigs, center_eigs, unit_eigs, stable_eigs, classify_eigs
export stability_index
export wrapto2pi, wraptoperiod, wraptopi, wrapto360, wrapto180
export to_dict, to_mat
export getpatchpoints

include("Orbits/manifolds.jl")
export subspace_stepoff, stable_manifold, unstable_manifold, linear_invariant_curve_2d

include("Orbits/quasiperiodic_orbit.jl")
export QuasiPeriodicOrbit
export reftraj,invariantcurve, DGmat, x2u, u2x

##################################################################
######################### CONSTRAINTS ############################
##################################################################

include("Targeting/Constraints/continuity_constraint.jl")
export ContinuityConstraint
export x1, x2, tof, dm
export evalconstraint, partials, cctspan, solvefunction, tangentsolvefunction

include("Targeting/Constraints/tofconstraint.jl")
export TOFConstraint
export td, tvec

include("Targeting/Constraints/invariance_constraint_2d.jl")
export InvarianceConstraint2D
export u0, u0vec, u0svec, xstar, strobetime, rotationangle, rhoval, kvec, thvec, Dmat, Qmat
export invariant_rotation, propagate_invariant_curve, ignore_imag

include("Targeting/Constraints/jacobi_constraint.jl")
export JacobiConstraint

include("Targeting/Constraints/q2p_rcont_constraint.jl")
export Q2P_PositionContinuityConstraint
export q2p_qpo, q2p_po, q2p_tht, q2p_thr, q2p_tau

include("Targeting/Constraints/perpendicular_crossing_constraint.jl")
export PerpendicularCrossingConstraint
export pcc_x1, pcc_T, pcc_tmax, xz_plane_crossing

include("Targeting/Constraints/continuity_epoch_constraint.jl")
export ContinuityEpochConstraint

include("Frames/reference_frames.jl")
export ReferenceFrame
export EM_BCR, EM_ECR, EM_MCR, EM_LVLH, EM_ICR, EM_TCR, EM_ECAI, EM_MCAI, EM_VCR, EM_VNC, EM_RECAI, EM_RMCAI, EJ2K, MJ2K
export isrelativeframe, isinertialframe
export frameconvert

include("Relative_Dynamics/lvlh.jl")
export LVLHModel

include("TwoBody/orbital_elements.jl")
export ClassicalOrbitalElements, sma, ecc, raan, inc, aop

include("io.jl")
export from_dict, from_mat

include("matlab_compat.jl")
export MatlabODESolver, matlabsolve

end
