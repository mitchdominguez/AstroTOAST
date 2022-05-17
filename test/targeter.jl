using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test


@testset "constraint.jl" begin
    # Define some free variables
    X1 = FreeVariable("x1",[0.72, 0, 0.71, 0, 0.18, 0])
    X2 = FreeVariable("x2",[0.72, 0, -0.71, 0, 0.18, 0])
    X3 = FreeVariable("x3",[0.72, 0, 0.71, 0, 0.18, 0])
    X4 = FreeVariable("x4",[0.72, 0, -0.71, 0, 0.18, 0])
    X5 = FreeVariable("x5",[0.72, 0, 0.71, 0, 0.18, 0])
    T = FreeVariable("T", 3.0, false) # inactive time
    # T = FreeVariable("T", 3.0, true) # active time
    # Tvec = Vector{FreeVariable}()
    # push!(Tvec,FreeVariable("T1", 3.0, true))# active time
    # push!(Tvec,FreeVariable("T2", 3.0, true))# active time
    # push!(Tvec,FreeVariable("T3", 3.0, true))# active time
    # push!(Tvec,FreeVariable("T4", 3.0, true))# active time

    # Define a free variable vector
    xv = XVector(X1, X2, X3, X4, X5, T)
    # xv = XVector(X1, X2, X3, X4, X5, Tvec...)

    # Define the model
    model = Cr3bpModel(Bodies["Earth"],Bodies["Moon"])

    # Define some constraints
    cc1 = ContinuityConstraint(X1, X2, T, model)
    cc2 = ContinuityConstraint(X2, X3, T, model)
    cc3 = ContinuityConstraint(X3, X4, T, model)
    cc4 = ContinuityConstraint(X4, X5, T, model)
    # cc1 = ContinuityConstraint(X1, X2, Tvec[1], model)
    # cc2 = ContinuityConstraint(X2, X3, Tvec[2], model)
    # cc3 = ContinuityConstraint(X3, X4, Tvec[3], model)
    # cc4 = ContinuityConstraint(X4, X5, Tvec[4], model)

    # Define a FXVector
    fx = FXVector(cc1, cc2, cc3, cc4)

    # Test initial norm
    @test norm(fx) == 0.2709052919898982

    # Create targeter object
    maxiter = 10
    tol = 1e-12
    targ = Targeter(xv, fx, maxiter, tol)

    # Test DFX matrix
    funcmat = [AstroTOAST.__dCC_dx1{6}()  AstroTOAST.__dCC_dx2{6}()  AstroTOAST.__NP{6}()       AstroTOAST.__NP{6}()       AstroTOAST.__NP{6}();
               AstroTOAST.__NP{6}()       AstroTOAST.__dCC_dx1{6}()  AstroTOAST.__dCC_dx2{6}()  AstroTOAST.__NP{6}()       AstroTOAST.__NP{6}();
               AstroTOAST.__NP{6}()       AstroTOAST.__NP{6}()       AstroTOAST.__dCC_dx1{6}()  AstroTOAST.__dCC_dx2{6}()  AstroTOAST.__NP{6}();
               AstroTOAST.__NP{6}()       AstroTOAST.__NP{6}()       AstroTOAST.__NP{6}()       AstroTOAST.__dCC_dx1{6}()  AstroTOAST.__dCC_dx2{6}()]

    @test DFX(targ) == funcmat


    # Test targeting
    errhist = [ 0.2709052919898982;
               0.03828309987199947;
               0.0012477184427681029;
               1.505021139162566e-6;
               1.9517145262023064e-12;
               1.1787495256637557e-13 ]
    
    Xhist, err = AstroTOAST.target(targ);

    @test err == errhist
    @test size(err) == size(Xhist)
    @test size(Xhist) == (6,)

    # Test norm(fx, xvorig, xveval)
    for i=1:length(Xhist)
        @test norm(fx, xv, Xhist[i]) == errhist[i]
    end
    @test tovector(Xhist[1][1]) == [0.72, 0, 0.71, 0, 0.18, 0]


end
