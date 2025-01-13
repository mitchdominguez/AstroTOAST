using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test
using MAT
using MATLAB
using SPICE

include("display.jl")
load_only_default_kernels()

@testset "test_ephemeris_rotation.jl" begin

    debug = false

    ls_ricardo = 3.8474799201129237e+05
    ts_ricardo = 3.756998590849912e+05
    mstar = dimensional_mass(AstroTOAST.em_cr3bp)

    dquant_atd = DimensionalQuantitySet(; mass=mstar, length=ls_ricardo, time=ts_ricardo)

    vars = matread("test/L2Halo9To2_5year_Jan22.mat")

    utcepoch = TDBEpoch(2023, 1, 1, 12, 0, 0.0)
    global hfem_atd = HFEModel(AstroTOAST.Moon, dquant_atd, [AstroTOAST.Earth, AstroTOAST.Sun], utcepoch)


    ind1 = 2
    ind2 = 3
    
    q0_atd = vec(vars["data"]["nodes"]["node"][ind1])
    q1_atd = vec(vars["data"]["nodes"]["node"][ind2])

    t0_nd = vars["data"]["t0s"]["t0"][ind1]
    t1_nd = vars["data"]["t0s"]["t0"][ind2]

    dt0_atd = vars["data"]["dts"]["dt"][ind1]

    if debug
        println("TOF (days): ", dt0_atd*dimensional_time(dquant_atd)*sec2day)

        # hfem_atd = HFEModel(AstroTOAST.Moon, dquant_atd, [AstroTOAST.Earth], utcepoch)


        println("Acceleration at t0: ", hfem_atd(q0_atd, t0_nd; debug=true), "\n")
        # println("Acceleration at tf: ", hfem_atd(q1_atd, t1_nd; debug=true), "\n")
    end

    sol = solve(hfem_atd, SVector{6,Float64}(q0_atd), (t0_nd, t1_nd))

    if debug

        println("Final state  (me): ", sol.u[end])
        println("Final state (ATD): ", q1_atd)
        println("Error (dimensional): ", dimensionalize_state(hfem_atd, sol.u[end] - q1_atd))
        println("Position Error (km): ", norm(dimensionalize_state(hfem_atd, sol.u[end] - q1_atd)[1:3]))
        println("Velocity Error (km/s): ", norm(dimensionalize_state(hfem_atd, sol.u[end] - q1_atd)[4:6]))

        println("\nNRHO Final Position (me):  ", frameconvert(sol.u[begin], to_ephemeris_time(t1_nd, hfem_atd), MJ2K(), EM_BCR()))
        println("NRHO Final Position (ATD): ", frameconvert(q0_atd, to_ephemeris_time(t1_nd, hfem_atd), MJ2K(), EM_BCR()))

        println("\nNondimensional state error: ", norm(sol.u[end]-q1_atd))
    end

    @test norm(sol.u[end]-q1_atd) < 1e-9


    # global nrhosol = sol

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################
    # plotstr = """
    # figure
    # hold on
    # grid on
    # axis equal
    # xlabel("x [ndim]")
    # ylabel("y [ndim]")
    # zlabel("z [ndim]")
    # plot3(1-0.012150584269940356,0,0, 'k.','markers',20)
    # title("EM Pulsating-Rotating")
    # """
    # eval_string(sesh, plotstr)


    # for i = 1:50
    # # for i = 2:2
        # println("i\n---")
        # q0_atd = vec(vars["data"]["nodes"]["node"][i])
        # q1_atd = vec(vars["data"]["nodes"]["node"][i+1])

        # t0_nd = vars["data"]["t0s"]["t0"][i]
        # t1_nd = vars["data"]["t0s"]["t0"][i+1]

        # sol = solve(hfem_atd, SVector{6,Float64}(q0_atd), (t0_nd, t1_nd))

        # hfemtraj = Trajectory(hfem_atd, sol.u[begin], [t0_nd, t1_nd])

        # # println(get_epoch(hfem_atd))

        # atd_proptimes = LinRange(tspan(hfemtraj)..., 300)

        # # println("-----   ", atd_proptimes[begin], "  ", atd_proptimes[end])

        # q_nrho_atd = hfemtraj(atd_proptimes)

        # R_q_nrho_atd = frameconvert(q_nrho_atd, to_ephemeris_time(atd_proptimes, hfem_atd), MJ2K(), EM_BCR())

        # # println(R_q_nrho_atd[1])

        # R_q_nrho_atd_mat = hcat(R_q_nrho_atd...)

        # put_variable(sesh, :R_qhist_nrho, R_q_nrho_atd_mat)

        # plotstr = """
        # plot3(R_qhist_nrho(1,:), R_qhist_nrho(2,:), R_qhist_nrho(3,:), 'k-','linewidth',2)
        # xlabel("x [ndim]")
        # ylabel("y [ndim]")
        # zlabel("z [ndim]")
        # title("EM Pulsating-Rotating")
        # """
        # eval_string(sesh, plotstr)
    # end

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################



    x_embcr = [0.884916533745800,0,0,0,0.470613569321843,0]
    utcepoch = UTCEpoch(2025, 1, 1, 12, 0, 0.0)

    model = AstroTOAST.em_cr3bp

    global hfem = HFEModel(AstroTOAST.Moon, dimensional_quantity_set(model), [AstroTOAST.Earth, AstroTOAST.Sun], utcepoch)

    I_R_ME_spice = ephemeris_state(;target=AstroTOAST.Earth, epoch=str2et(utcepoch), frame="J2000", observer=AstroTOAST.Moon)

    x_mj2k = frameconvert(x_embcr, utcepoch, EM_BCR(), MJ2K())

    # x_mj2k_nadia = [-0.051220074685291,0.077314030343455,0.041924988067173,0.320729606562720,0.160572871492497,0.086751674764737]
    #
    global x_mj2k_nadia = [-0.0511737477981375
                    0.0772441023243908
                    0.0418870682827328
                    0.320874752694711
                    0.160645538719728
                    0.0867909342212326]

    # println(x_embcr)
    # println(x_mj2k)
    # println(x_mj2k_nadia)
    # println(x_mj2k_nadia - x_mj2k)
    # println(I_R_ME_spice)
    
    acc_me = hfem(x_mj2k_nadia, 0.0; debug=false)
    acc_nadia = [0.320874752694711, 0.160645538719728, 0.0867909342212326, 0.467751670849167, -0.705001755076245, -0.382467436423266]

    # println("Acceleration (Mitch): ", acc_me)
    # println("Acceleration (Nadia): ", acc_nadia)
    # println("Acceleration Difference: ", acc_nadia - acc_me)

    @test norm(acc_nadia - acc_me) <= 1e-12

    sol_DRO = solve(hfem, SVector{6,Float64}(x_mj2k_nadia), (0, 6.908*3))

    global dro_mj2k = sol_DRO.u

    tsol_DRO = tangent_solve(hfem, SVector{6,Float64}(x_mj2k_nadia), (0, 6.908*3))

    global dro_stm_mj2k = tsol_DRO.u


    @test norm(dro_stm_mj2k[end][:,end] - [  0.1425614207307463
                                           -0.05961287646580876
                                           -0.03153884264446574
                                           0.044483667864847386
                                           -0.33391360838597833
                                           -0.18837817987835812
                                          ]) <= 1e-11


    # @test norm(dro_stm_mj2k[end][:,end] - [0.14256142072712272
                                           # -0.0596128764642774
                                           # -0.03153884264364764
                                           # 0.04448366786858148
                                           # -0.333913608376461
                                           # -0.1883781798729954]) <= 1e-12
end
