using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test
using BenchmarkTools

######## Obtain QPO and PO ########
include("nrho92.jl");
include("quasiperiodic_orbit.jl");

######## Define FreeVariables ########
# thT = FreeVariable("thT", 0., includeinds=[])
thT = FreeVariable("thT", 0.001710)
thR = FreeVariable("thR", 0.363052)
tau = FreeVariable("τ", 3.114937396)
xv = XVector(thT, thR, tau)

######## Define Constraint ########
q2p = Q2P_PositionContinuityConstraint(qpo92, nrho92, thT, thR, tau)
fx = FXVector(q2p)

######## Target ########
maxiter = 100
tol = 1e-12
targ = Targeter(xv, fx, maxiter, tol);
Xhist, err = target(targ,debug=true)
err

