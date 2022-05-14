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
export dimension, name, value, update!, update, active
export XVector
export numels, tovector, tosvector

include("Targeting/constraint.jl")
include("Targeting/continuity_constraint.jl")
export ContinuityConstraint
export x1, x2, tof, dm
export evalconstraint, partials, cctspan

end
