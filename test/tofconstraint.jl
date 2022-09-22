using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test

# Define the model
model = Cr3bpModel(Bodies["Earth"],Bodies["Moon"])

# Define some free variables
# X1 = FreeVariable("x1", [0.987, 0, 0.008, 0, 1.667, 0], [2,4,6]) # fix x0, y0, xdot0, zdot0
X1 = FreeVariable("x1", [0.987, 0, 0.008, 0, 1.667, 0]) # Don't fix any variables
X2 = FreeVariable("x2", [0.988, 0.023, -0.009,0.079, 0.488, -0.802])
X3 = FreeVariable("x3", [1.021, 0.011, -0.178,0.015, -0.100, -0.058])


T1 = FreeVariable("T1", 0.0239)
T2 = FreeVariable("T2", 0.6156)
T3 = FreeVariable("T3", 0.8718)

# Define XVector
xv = XVector(X1, X2, X3, T1, T2, T3)

# Define constraints

#   keep full state in X, remove some constraints
cc = Vector{ContinuityConstraint}()
push!(cc, ContinuityConstraint(X1, X2, T1, model))
push!(cc, ContinuityConstraint(X2, X3, T2, model))
push!(cc, ContinuityConstraint(X3, X1, T3, model, includeinds=[1,2,3,4,6])) # Remove ydot constraint

# Define FXVector
fx = FXVector(cc...) # FX vector for full X, rm in FX

# Define Targeter
maxiter = 20
tol = 1e-12
targ = Targeter(xv, fx, maxiter, tol);

# Target fx rc version
Xhist, err = target(targ);


@testset "tofconstraint.jl" begin
    Td = 1.511261560928471 # 9:2 NRHO
    # Define time of flight constraint
    tofc = TOFConstraint(Td, T1, T2, T3)

    # Evaluate the TOF constraint
    @test td(tofc) == Td
    @test tvec(tofc) == [T1, T2, T3]
    @test evalconstraint(tofc)≈[0.023002456816265004] atol=1e-10

    # Make a new FXVector to add in the TOF constraint
    fx2 = FXVector(cc..., tofc)
    @test numels(fx2) == 4
    @test fx2[4] == tofc
    
    # Reset XVector
    AstroTOAST.update!(xv, Xhist[1])

    # Target with tof constraint
    targ2 = Targeter(xv, fx2, maxiter, tol)
    Xhist2, err2 = target(targ2)

    # Results
    @test length(err2) == 8
    @test norm(tofullvector(fx2)) ≈ 6.459955702918004e-14 atol=1e-12
    @test evalconstraint(tofc)[1] ≈ 0.0 atol=1e-14


end
