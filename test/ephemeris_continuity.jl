using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test
using MAT
using MATLAB
using SPICE
using SparseArrays

include("display.jl")
load_only_default_kernels()

debug = false

@testset "test_ephemeris_singleshooter.jl" begin

    debug = false

    ## Specify epoch 
    utcepoch = UTCEpoch(2025, 1, 1, 12, 0, 0.0)

    model = AstroTOAST.em_cr3bp
    lstar = dimensional_length(model)
    tstar = dimensional_time(model)

    ## 1 day, in nondimensional time
    day_ndim = day2sec/tstar

    global hfem = HFEModel(AstroTOAST.Moon, dimensional_quantity_set(model), [AstroTOAST.Earth, AstroTOAST.Sun], utcepoch)

    global x_mj2k_nadia = [-0.0511737477981375, 0.0772441023243908, 0.0418870682827328, 0.320874752694711, 0.160645538719728, 0.0867909342212326]
    sol_DRO = solve(hfem, SVector{6,Float64}(x_mj2k_nadia), (0, 9day_ndim))

    ## Patch points
    X1 = FreeVariable("X1", sol_DRO(0day_ndim), includeinds=1:6)
    X2 = FreeVariable("X2", sol_DRO(3day_ndim).+[0.0005,0,0,0,0.0001,0], includeinds=1:6)
    X3 = FreeVariable("X3", sol_DRO(6day_ndim).+[0.0,0.00001,0,0,0,-0.0002], includeinds=1:6)
    X4 = FreeVariable("X4", sol_DRO(9day_ndim).+[0.0005,0,0,0,0,0], includeinds=1:6)

    ## Epochs corresponding to each patch point
    T1 = FreeVariable("T1", 0day_ndim, includeinds=1:1)
    T2 = FreeVariable("T2", 3.01day_ndim, includeinds=1:1)
    T3 = FreeVariable("T3", 6day_ndim, includeinds=1:1)
    T4 = FreeVariable("T4", 9day_ndim, includeinds=1:1)

    ## Times of flight
    TOF1 = FreeVariable("TOF1", 3day_ndim, includeinds=1:1)
    TOF2 = FreeVariable("TOF2", 3day_ndim, includeinds=1:1)
    TOF3 = FreeVariable("TOF3", 3day_ndim, includeinds=1:1)
    TOF4 = FreeVariable("TOF4", 3day_ndim, includeinds=1:1)

    ###### Try a simple single shooter

    ## XVector
    xv = XVector(X1, TOF1, T1, T2)

    ## Constraints
    cc = ContinuityEpochConstraint(X1, X2, T1, T2, TOF1, hfem)
    fx = FXVector(cc)

    targ = Targeter(xv, fx, 20, 1e-12)


    if debug
        println("Initial state discontinuity: ", ignore_zero(evalconstraint(cc)[1:6]))
        println("Initial time discontinuity: ", T1[1] + TOF1[1] - T2[1])
    end

    Xhist, err = target(targ; debug=false)
    xhistmat = hcat(tofullvector.(Xhist)...)

    if debug
        println("Final state discontinuity: ", ignore_zero(evalconstraint(cc)[1:6]))
        println("Final time discontinuity: ", evalconstraint(cc)[7])
        println("Final time discontinuity: ", T1[1] + TOF1[1] - T2[1])
    end

    @test norm(evalconstraint(cc)[1:6]) < 1e-12
    @test norm(evalconstraint(cc)[7]) < 1e-12

end

@testset "test_ephemeris_multipleshooter.jl" begin

    ## Specify epoch 
    utcepoch = UTCEpoch(2025, 1, 1, 12, 0, 0.0)

    model = AstroTOAST.em_cr3bp
    lstar = dimensional_length(model)
    tstar = dimensional_length(model)

    ## 1 day, in nondimensional time
    day_ndim = day2sec/tstar

    global hfem = HFEModel(AstroTOAST.Moon, dimensional_quantity_set(model), [AstroTOAST.Earth, AstroTOAST.Sun], utcepoch)

    global x_mj2k_nadia = [-0.0511737477981375, 0.0772441023243908, 0.0418870682827328, 0.320874752694711, 0.160645538719728, 0.0867909342212326]
    sol_DRO = solve(hfem, SVector{6,Float64}(x_mj2k_nadia), (0, 9day_ndim))

    ## Patch points
    X1 = FreeVariable("X1", sol_DRO(0day_ndim), includeinds=1:6)
    X2 = FreeVariable("X2", sol_DRO(3day_ndim).+[0.0005,0,0,0,0.0001,0], includeinds=1:6)
    X3 = FreeVariable("X3", sol_DRO(6day_ndim).+[0.0,0.00001,0,0,0,-0.0002], includeinds=1:6)
    X4 = FreeVariable("X4", sol_DRO(9day_ndim).+[0.0005,0,0,0,0,0], includeinds=1:6)

    ## Epochs corresponding to each patch point
    T1 = FreeVariable("T1", 0day_ndim, includeinds=1:1)
    T2 = FreeVariable("T2", 3.01day_ndim, includeinds=1:1)
    T3 = FreeVariable("T3", 6day_ndim, includeinds=1:1)
    T4 = FreeVariable("T4", 9day_ndim, includeinds=1:1)

    ## Times of flight
    TOF1 = FreeVariable("TOF1", 3day_ndim, includeinds=1:1)
    TOF2 = FreeVariable("TOF2", 3day_ndim, includeinds=1:1)
    TOF3 = FreeVariable("TOF3", 3day_ndim, includeinds=1:1)

    ###### Try a multiple shooter

    ## XVector
    # xv = XVector(X1, X2, TOF1, TOF2, T1, T2)
    # xv = XVector(X1, X2, X3, TOF1, TOF2, TOF3, T1, T2, T3)
    xv = XVector(X1, X2, X3, X4, TOF1, TOF2, TOF3, T1, T2, T3, T4)

    ## Constraints
    cc1 = ContinuityEpochConstraint(X1, X2, T1, T2, TOF1, hfem)
    cc2 = ContinuityEpochConstraint(X2, X3, T2, T3, TOF2, hfem)
    cc3 = ContinuityEpochConstraint(X3, X4, T3, T4, TOF3, hfem)
    fx = FXVector(cc1, cc2, cc3)

    global targ = Targeter(xv, fx, 20, 1e-12)

    if debug
        println("Initial state discontinuities: ", [norm(evalconstraint(x)[1:6]) for x in [cc1,cc2,cc3]])
        println("Initial time discontinuities: ", [evalconstraint(x)[end] for x in [cc1,cc2,cc3]])
    end

    __sparsedf__(T::Targeter) = sparse(evalDFXMatrix(T))

    # Xhist, err = target(targ; debug=false)
    # @time Xhist, err = target(targ; debug=false, inversion_method=:fancy)
    # @time Xhist, err = target(targ; debug=false, inversion_method=:fancy, eval_DF_func=__sparsedf__)
    # @time Xhist, err = target(targ; debug=false, inversion_method=:fancy)
    @time Xhist, err = target(targ; debug=false, inversion_method=:backslash)
    xhistmat = hcat(tofullvector.(Xhist)...)

    if debug
        println("Initial state discontinuities: ", [norm(evalconstraint(x)[1:6]) for x in [cc1,cc2,cc3]])
        println("Initial time discontinuities: ", [evalconstraint(x)[end] for x in [cc1,cc2,cc3]])
    end

    for cc in [cc1, cc2, cc3]
        @test norm(evalconstraint(cc)[1:6]) < 1e-12
        @test norm(evalconstraint(cc)[7]) < 1e-12
    end

    ephem_traj = Trajectory(hfem, [X1, X2, X3], [T1, T2, T3, T4])
    @test tspan(ephem_traj)[1] == T1[1]
    @test tspan(ephem_traj)[2] == T4[1]

end
