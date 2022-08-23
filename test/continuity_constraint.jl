using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test


@testset "continuity_constraint.jl" begin
    # Define some free variables
    X1 = FreeVariable("x1",[0.9, 0, 0, 0, -0.4, 0])
    X2 = FreeVariable("x2", [0.95, 0.05, 0, 0.48, -0.2, 0])
    X3 = FreeVariable("x2", [0.85, 0.85, 3, 1.48, -0.9, -2])
    T = FreeVariable("T", 2.0)
    Tf = FreeVariable("T", 2.0, includeinds=[]) # inactive time

    # Define a free variable vector
    xv = XVector(X1, X2, T)
    Xvf = XVector(X1, X2, Tf)

    # Define the model
    model = Cr3bpModel(Bodies["Earth"],Bodies["Moon"])

    # constructor
    cc = ContinuityConstraint(X1, X2, T, model)
    ccf = ContinuityConstraint(X1, X2, Tf, model)
    @test_throws MethodError ContinuityConstraint(X1, X2, 2.0, model)

    # x1, x2, tof, model
    @test x1(cc) == X1
    @test x2(cc) == X2
    @test tof(cc) == T
    @test dm(cc) == model
    @test x1(ccf) == X1
    @test x2(ccf) == X2
    @test tof(ccf) == Tf
    @test dm(ccf) == model

    # cctspan
    @test cctspan(cc) == (0, 2.0)
    @test cctspan(ccf) == (0, 2.0)

    # evalconstraint
    constval = @SVector [0.005382188230026896, 0.0024307725939954616, 0.0, 
                         -0.0008239789557616395, 0.03766157082338922, 0.0]
    @test evalconstraint(cc)≈constval atol=1e-12


    # partials
    @test partials(cc, X1) == AstroTOAST.__dCC_dx1{6}()
    @test partials(cc, X2) == AstroTOAST.__dCC_dx2{6}()
    @test partials(cc, X3) == AstroTOAST.__NP{6}()
    @test partials(cc, T) == AstroTOAST.__dCC_dt{1}()
    @test partials(ccf, Tf) == AstroTOAST.__NP{1}()


    # variables in an XVector still result in the right partial function being called
    @test partials(cc, xv[1]) == AstroTOAST.__dCC_dx1{6}()
    @test partials(cc, xv[2]) == AstroTOAST.__dCC_dx2{6}()
    @test partials(cc, xv[3]) == AstroTOAST.__dCC_dt{1}()

    # partial derivative functions
    @test AstroTOAST.__dCC_dx1{6}()(cc) ≈ 
    [187.566    22.6487  0.0         19.5816    32.929    0.0
     -113.84    -14.2254  0.0        -12.0144   -19.9858   0.0
     0.0       0.0     0.153128     0.0        0.0     -0.189341
     833.483   104.65    0.0         88.5242   145.504    0.0
     -1273.21   -150.901   0.0       -131.465   -224.878    0.0
     0.0       0.0     5.13867      0.0        0.0      0.17656] atol=1e-2

    @test AstroTOAST.__dCC_dx2{6}()(cc) == -I(6)
    @test AstroTOAST.__NP{6}()(cc) == zeros(6,6)

    constvec = [0.4791760210447607, -0.16233842917687474, 0.0, 
                1.2621101313341967, -3.679141924348646, -0.0] 
    @test AstroTOAST.__dCC_dt{1}()(cc) ≈ constvec atol=1e-10

end
