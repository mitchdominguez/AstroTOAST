### Test propagation and targeting in the LVLH frame
using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test

# Use the 9:2 NRHO
include("nrho92.jl")

# Initialize model
model = Cr3bpModel(Bodies["Earth"],Bodies["Moon"]) 
Lmod = LVLHModel(model)
# nrho92 = generate_nrho92(model,return_targeter=false)

# Initial conditions
xt = nrho92(pi)
xc = nrho92(pi+0.01)
q0 = vcat(frameconvert(xt, xc, EM_BCR(), EM_LVLH())...)

sol_init = solve(Lmod, q0_conv, period(nrho92); abstol=AstroTOAST.DEFAULT_ABS_TOL)

@testset "lvlh_propagation.jl" begin
    @test norm(sol_init[end]-q0 ) <= AstroTOAST.DEFAULT_CONVERGENCE_TOL
end

# Set up targeting problem
X1 = FreeVariable("X1", q0; includeinds=10:12)
X2 = FreeVariable("X2", [xt..., zeros(6)...]; includeinds=10:12)
T = FreeVariable("tof", 8*hr2sec/dimensional_time(Lmod); includeinds=[])
xv = XVector(X1, X2, T)
cc = ContinuityConstraint(X1, X2, T, Lmod;includeinds=7:12)
fx = FXVector(cc)
# targ = Targeter(xv, fx, maxiter, tol)

# Xhist, err = target(targ)
q0_conv = tofullvector(X1)

sol = solve(Lmod, q0_conv, cctspan(cc); abstol=AstroTOAST.DEFAULT_ABS_TOL)

println(sol[end][7:12])

@testset "lvlh_targeting.jl" begin

    @test err[end] < AstroTOAST.DEFAULT_CONVERGENCE_TOL
    @test norm(sol[end][7:9]) <= AstroTOAST.DEFAULT_CONVERGENCE_TOL

    

end
