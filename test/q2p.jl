using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test
using BenchmarkTools

######## Obtain QPO and PO ########
include("nrho92.jl");
include("quasiperiodic_orbit.jl");

######## Define FreeVariables ########
thT = FreeVariable("thT", 0., includeinds=[])
thR = FreeVariable("thR", 0.)
tau = FreeVariable("τ", Float64(pi))
xv = XVector(thT, thR, tau)

######## Define Constraint ########
q2p = Q2P_PositionContinuityConstraint(qpo92, nrho92, thT, thR, tau)
fx = FXVector(q2p)

######## Target ########
maxiter = 100
tol = 1e-12
targ = Targeter(xv, fx, maxiter, tol);

# Target Position Correspondence
ths = LinRange(0, 2pi, 151)
for th in ths
    update!(thT, [th])
    update!(tau, [th])
    try
        target(targ,debug=false);
        println("$(th) DID WORK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    catch
        println("$(th) did not work")
    end
end
