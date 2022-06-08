using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test


@testset "poms.jl" begin
    # TARGET ORBITS NEAR THE 9:2 NRHO

    # Define the model
    model = Cr3bpModel(Bodies["Earth"],Bodies["Moon"])

    # Define some free variables
    # X1 = FreeVariable("x1", [0.987, 0, 0.008, 0, 1.667, 0], [2,4,6]) # fix x0, y0, xdot0, zdot0
    X1 = FreeVariable("x1", [0.987, 0, 0.008, 0, 1.667, 0], includeinds=2:6) # fix y0
    X2 = FreeVariable("x2", [0.988, 0.023, -0.009,0.079, 0.488, -0.802])
    X3 = FreeVariable("x3", [1.021, 0.011, -0.178,0.015, -0.100, -0.058])

    X1_full = FreeVariable("x1", [0.987, 0, 0.008, 0, 1.667, 0]) # Don't fix any variables

    T1 = FreeVariable("T1", 0.0239)
    T2 = FreeVariable("T2", 0.6156)
    T3 = FreeVariable("T3", 0.8718)

    # Define XVector
    xv = XVector(X1, X2, X3, T1, T2, T3)
    xv_full = XVector(X1_full, X2, X3, T1, T2, T3)

    # Define constraints
    #   removing states in X, keeping all constraints
    cc_rxfc = Vector{ContinuityConstraint}()
    push!(cc_rxfc, ContinuityConstraint(X1, X2, T1, model))
    push!(cc_rxfc, ContinuityConstraint(X2, X3, T2, model))
    push!(cc_rxfc, ContinuityConstraint(X3, X1, T3, model))

    #   keep full state in X, remove some constraints
    cc_fxrc = Vector{ContinuityConstraint}()
    push!(cc_fxrc, ContinuityConstraint(X1_full, X2, T1, model))
    push!(cc_fxrc, ContinuityConstraint(X2, X3, T2, model))
    push!(cc_fxrc, ContinuityConstraint(X3, X1_full, T3, model, includeinds=[1,2,3,4,6])) # Remove ydot constraint

    # Define FXVector
    fx_rxfc = FXVector(cc_rxfc...) # FX vector for rm in X, full FX
    fx_fxrc = FXVector(cc_fxrc...) # FX vector for full X, rm in FX

    # Define Targeter
    maxiter = 20
    tol = 1e-12
    targ_fxrc = Targeter(xv_full, fx_fxrc, maxiter, tol);

    # Target fx rc version
    Xhist, err = target(targ_fxrc);
    
    # Test that norm of vector and fullvector are different, but that both satisfy the tolerance
    @test norm(fx_fxrc) == err[end]
    @test norm(tovector(fx_fxrc)) == err[end]
    @test norm(tofullvector(fx_fxrc)) != norm(fx_fxrc)
    @test norm(fx_fxrc) ≈ 6.72863818142289e-14 atol=AstroTOAST.DEFAULT_CONVERGENCE_TOL
    @test norm(tofullvector(fx_fxrc)) ≈ 7.521894782093819e-14 atol=AstroTOAST.DEFAULT_CONVERGENCE_TOL
    @test norm(tofullvector(fx_fxrc)) < tol

    # Initial error
    @test err[1] ≈ 0.7464423331586798 atol=AstroTOAST.DEFAULT_CONVERGENCE_TOL
    @test norm(fx_fxrc, xv_full, Xhist[1]) ≈ 0.7464423331586798 atol=AstroTOAST.DEFAULT_CONVERGENCE_TOL# Initial error

    # Test number of iterations (this was a bad initial guess)
    @test length(err) == 9

    # Reset the test to run the alternate example where a free variable is fixed, and we use all the constraints
    update!(xv, Xhist[1])

    @test norm(fx_rxfc) ≈ 0.8417710513141402 atol=AstroTOAST.DEFAULT_CONVERGENCE_TOL

    targ_rxfc = Targeter(xv, fx_rxfc, maxiter, tol);
    Xhist, err = target(targ_rxfc);

    # Number of iterations
    @test length(err) == 15

    # Test that norm of vector and fullvector are the same, since no constraints are removed
    @test norm(fx_rxfc) == err[end]
    @test norm(tovector(fx_rxfc)) == err[end]
    @test norm(tofullvector(fx_rxfc)) == norm(fx_rxfc)
    @test norm(tofullvector(fx_rxfc)) < tol


end
