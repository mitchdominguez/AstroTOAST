using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test


@testset "poms.jl" begin
    # TARGET AN ORBIT NEAR THE 9:2 NRHO

    # Define the model
    model = Cr3bpModel(Bodies["Earth"],Bodies["Moon"])

    # Define some free variables
    X1 = FreeVariable("x1", [0.987, 0, 0.008, 0, 1.667, 0], [2,4,6]) # fix x0, y0, xdot0, zdot0
    # X1 = FreeVariable("x1", [0.987, 0, 0.008, 0, 1.667, 0], [2]) # fix y0
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
    push!(cc_rxfc, ContinuityConstraint(X1, X2, T1, model, [2,4,6]))
    push!(cc_rxfc, ContinuityConstraint(X2, X3, T2, model))
    push!(cc_rxfc, ContinuityConstraint(X3, X1, T3, model))

    #   keep full state in X, remove some constraints
    cc_fxrc = Vector{ContinuityConstraint}()
    push!(cc_fxrc, ContinuityConstraint(X1_full, X2, T1, model))
    push!(cc_fxrc, ContinuityConstraint(X2, X3, T2, model))
    push!(cc_fxrc, ContinuityConstraint(X3, X1_full, T3, model, 5)) # Remove ydot constraint

    # Define FXVector
    fx_rxfc = FXVector(cc_rxfc...) # FX vector for rm in X, full FX
    fx_fxrc = FXVector(cc_fxrc...) # FX vector for full X, rm in FX

    # Define Targeter
    maxiter = 15
    tol = 1e-12
    targ_fxrc = Targeter(xv_full, fx_fxrc, maxiter, tol);
    Xhist, err = target(targ_fxrc);
    
    println(err)

    # update!(xv, Xhist[end-4])
    update!(xv, Xhist[1])

    println(norm(fx_rxfc))

    targ_rxfc = Targeter(xv, fx_rxfc, maxiter, tol);
    Xhist, err = target(targ_rxfc);

    println(err)


end
