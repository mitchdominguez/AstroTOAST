using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test


@testset "constraint.jl" begin
    # Define some free variables
    X1 = FreeVariable("x1",[0.9, 0, 0, 0, -0.4, 0])
    X2 = FreeVariable("x2", [0.95, 0.05, 0, 0.48, -0.2, 0])
    X3 = FreeVariable("x2", [0.85, 0.85, 3, 1.48, -0.9, -2])
    T = FreeVariable("T", 2.0)
    Tf = FreeVariable("T", 2.0, false) # inactive time

    # Define a free variable vector
    xv = XVector(X1, X2, T)
    Xvf = XVector(X1, X2, Tf)

    # Define the model
    model = Cr3bpModel(Bodies["Earth"],Bodies["Moon"])

    # Define some constraints
    cc = ContinuityConstraint(X1, X2, T, model)
    ccf = ContinuityConstraint(X1, X2, Tf, model)
    cc2 = ContinuityConstraint(X2, X3, T, model)
    @test_throws MethodError ContinuityConstraint(X1, X2, 2.0, model)

    # Define a FXVector
    fx = FXVector(cc, cc2)

    # numels
    @test numels(fx) == 2

    # length
    @test length(fx) == 12

    # iterate
    @test iterate(fx,1) == (cc,2)
    @test iterate(fx,2) == (cc2,3)

    # tovector
    fx_evaled = [0.005382188230026896, 0.0024307725939954616, 0.0, 
                 -0.0008239789557616395, 0.03766157082338922, 0.0, 
                 0.04885903431274763, -0.8443322399829405, -3.0, 
                 -1.5427004376247462, 0.4847234176486228, 2.0]

    @test tovector(fx) == fx_evaled

    # tosvector
    @test tosvector(fx) == SVector(Tuple(fx_evaled))

    # getindex, x1, x2, tof, model
    @test x1(fx[1]) == X1
    @test x2(fx[1]) == X2
    @test tof(fx[1]) == T
    @test dm(fx[1]) == model
end
