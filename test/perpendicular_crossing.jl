using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test
using BenchmarkTools

### Target an L1 Lyapunov

# Dynamical model
model = Cr3bpModel(Bodies["Earth"],Bodies["Moon"])

# Obtain initial guess from the LVE at L₁
xi_0 = 0.01
index = 1
q0 = lyapunov_linear_ics(model, xi_0, 0, index) # Get linear stepoff from L1

############### FREEVARIABLE SETUP ############### 
# X1 = FreeVariable("x1", q0, includeinds=[1,3,5])
X1 = FreeVariable("x1", q0, includeinds=[5])
# X1 = FreeVariable("x1", q0)
T = FreeVariable("t", Float64(2π))

xv = XVector(X1, T)

############### CONSTRAINT SETUP ############### 
# pcc = PerpendicularCrossingConstraint(model, X1, T, includeinds=[1,2,4,5])
pcc = PerpendicularCrossingConstraint(model, X1, T, includeinds=[2,4])
fx = FXVector(pcc)

############### TARGETER SETUP ############### 
maxiter = 25
tol = 1e-12
targ = Targeter(xv, fx, maxiter, tol);

# Target perpendicular crossing
Xhist, err = target(targ,debug=false);

err

# @testset "perpendicular_crossing.jl" begin
    # @test 1==1
# end

