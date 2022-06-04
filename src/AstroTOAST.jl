module AstroTOAST

using StaticArrays
using ForwardDiff
using OrdinaryDiffEq
using DiffEqBase
using LinearAlgebra

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

include("Body/naifbody.jl")
export Bodies

include("CR3BP/cr3bp.jl")
export Cr3bpModel
export mass_ratio, primary_bodies
export distance_to_primary, primary_state
export pseudopotential, pseudopotential_gradient, pseudopotential_jacobian
export jacobi_constant

include("Targeting/free_variable.jl")
export FreeVariable
export full_length, name, fullvalue, tofullvector, value, update!, update, active, removeinds
export XVector
export numels, tovector, tosvector, tofullsvector

include("Targeting/constraint.jl")
export FXVector

include("Targeting/targeter.jl")
export Targeter
export X, FX, DFX, getmaxiter, gettol, evalDFXMatrix, target

include("Targeting/Constraints/continuity_constraint.jl")
export ContinuityConstraint
export x1, x2, tof, dm
export evalconstraint, partials, cctspan

include("Targeting/Constraints/tofconstraint.jl")
export TOFConstraint
export td, tvec

include("Orbits/trajectory.jl")
export Trajectory
export x0, tspan, solvec, isperiodic, iscontinuous, stm

include("Orbits/periodic_orbit.jl")
export PeriodicOrbit
export family
export period, traj, monodromy, lambda, paireigs!, sorteigs!
export angle2time, time2angle
export unstable_eigs, center_eigs, unit_eigs, stable_eigs, classify_eigs
export stability_index

include("Orbits/manifolds.jl")
export subspace_stepoff, stable_manifold, unstable_manifold, linear_invariant_curve_2d

end
