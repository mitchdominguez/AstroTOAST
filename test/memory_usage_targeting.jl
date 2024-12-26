
using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test
using MAT
using MATLAB
using SPICE
using SparseArrays

include("display.jl")
load_only_default_kernels()

include("../../Formation_Flying/src/Orbits/ephemeris_transition.jl")
include("../../Formation_Flying/src/Constraints/continuity_epoch_T0_constraint.jl")

model = AstroTOAST.em_cr3bp
lstar = dimensional_length(model)
tstar = dimensional_time(model)

utcep = UTCEpoch(2025, 1, 1, 12, 0, 0.0)

halo21 = from_mat("../Formation_Flying/data/halo21_refnew.mat")[1]
nrho92 = from_mat("../Formation_Flying/data/reference_orbits.mat")[1]

x = halo21(0)

# disc_traj, longtarg = transition_to_ephemeris(halo21, utcep,
                                              # 600*period(halo21)*tstar*sec2day;
                                              # only_propagate_patch_points=true,
                                              # infinity_norm_tol=1e-12,
                                              # pts_per_rev=4,
                                              # variable_time=true,
                                              # tol=2.3e-11)


function get_nrho92_baseline(nrho92, utcep, tol, inftol, numrevs)
    disc_traj, longtarg = transition_to_ephemeris(nrho92, utcep,
                                                  numrevs*period(nrho92)*tstar*sec2day;
                                                  only_propagate_patch_points=false,
                                                  infinity_norm_tol=inftol,
                                                  pts_per_rev=3,
                                                  variable_time=true,
                                                  tol=tol)

    to_mat(disc_traj, "data/ephem_nrho92_$(numrevs)_revs")
end

tol2 = 2.5e-11
inftol = tol2*10
numrevs = 1200
get_nrho92_baseline(nrho92, utcep, tol2, inftol, numrevs)

# function test_df(longtarg)

# # evalDFXMatrix(longtarg)
# # @time evalDFXMatrix(longtarg)

# AstroTOAST.construct_sparse_DF(longtarg)
# @time AstroTOAST.construct_sparse_DF(longtarg)

# end

# test_df(longtarg)

# Xhist, err = target(longtarg; inversion_method=:fancy, infinity_norm_tol=1e-12, debug=true)
