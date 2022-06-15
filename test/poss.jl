using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test

# TARGET 2:1 Halo

# Define the model
model = Cr3bpModel(Bodies["Earth"],Bodies["Moon"])

# Define some free variables
X1 = FreeVariable("x1", [1.178977343109523, 0, -0.042926770700000, 0, -0.165685370741522, 0], includeinds=[1,3,5]) # 2:1 Halo
T1 = FreeVariable("T1", 3.400338564669278, includeinds=[])

# Define XVector
xv = XVector(X1, T1)

# Define constraints
cc = ContinuityConstraint(X1,X1,T1,model,includeinds=[2,4,6]) # constrain s.t we have a perpendicular crossing, must constrain y, xdot, zdot

# Define FXVector
fx = FXVector(cc) # FX vector for rm in X, full FX

# Define Targeter
maxiter = 20
tol = 1e-12
targ = Targeter(xv, fx, maxiter, tol);

# Target 2:1 halo
Xhist, err = target(targ);

@testset "poss.jl" begin

    @test norm(fx, xv, Xhist[1]) ≈ 0.0028364554191276227 atol=AstroTOAST.DEFAULT_CONVERGENCE_TOL
    @test norm(fx) ≈ 2.402323662518765e-13 atol=AstroTOAST.DEFAULT_CONVERGENCE_TOL
    @test norm(fx) < tol
    @test norm(fx) == err[end] < tol
    @test norm(tofullvector(fx)) ≈ 8.242143522346135e-13 atol=AstroTOAST.DEFAULT_CONVERGENCE_TOL

    @test length(err) == 5
end

## Test an overconstrained example. In this case, I am fixing both x0 and period, which overconstrains the system

# Define some free variables
X1_oc = FreeVariable("x1", [1.178977343109523, 0, -0.042926770700000, 0, -0.165685370741522, 0], includeinds=[3,5]) # 2:1 Halo

# Define XVector
xv_oc = XVector(X1_oc, T1)

# Define constraints
cc_oc = ContinuityConstraint(X1_oc,X1_oc,T1,model,includeinds=[2,4,6]) # constrain s.t we have a perpendicular crossing, must constrain y, xdot, zdot

# Define FXVector
fx_oc = FXVector(cc_oc) # FX vector for rm in X, full FX

# Target 
targ_oc = Targeter(xv_oc, fx_oc, maxiter, tol);

@testset "poss_overconstrained.jl" begin

    # Targeter fails
    @test_throws Exception Xhist_oc, err_oc = target(targ_oc, debug=false);
    @test_throws UndefVarError Xhist_oc
    @test_throws UndefVarError err_oc

end
