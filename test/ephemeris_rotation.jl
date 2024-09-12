using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test
using MAT
using SPICE

include("display.jl")
load_only_default_kernels()

@testset "test_ephemeris_rotation.jl" begin

    lstar = dimensional_length(AstroTOAST.em_cr3bp)
    tstar = dimensional_time(AstroTOAST.em_cr3bp)
    mu = AstroTOAST.em_cr3bp.μ
    GMstar = AstroTOAST.Moon.gravitational_parameter + AstroTOAST.Earth.gravitational_parameter
    epoch_tdb = TDBEpoch(2017, 8, 21, 0, 0, 0.0)

    ## Relevant vectors
    R_r_BE = [-mu,0,0,0,0,0]
    I_R_ME_spice = ephemeris_state(;target=AstroTOAST.Earth, epoch=str2et(epoch_tdb), frame="J2000", observer=AstroTOAST.Moon)
    I_r_ME_spice = nondimensionalize_state(AstroTOAST.em_cr3bp, I_R_ME_spice)

    R_r_BM = [1-mu,0,0,0,0,0]
    I_R_EM_spice = ephemeris_state(;target=AstroTOAST.Moon, epoch=str2et(epoch_tdb), frame="J2000", observer=AstroTOAST.Earth)
    I_r_EM_spice = nondimensionalize_state(AstroTOAST.em_cr3bp, I_R_EM_spice)

    ###### Test Rotating to MOONJ2000
    I_r_ME = frameconvert(R_r_BE, epoch_tdb, EM_BCR(), MJ2K())
    I_R_ME = dimensionalize_state(AstroTOAST.em_cr3bp, I_r_ME)
    println("Rotating -> MJ2K: ", norm(I_R_ME - I_R_ME_spice))
    @test norm(I_R_ME - I_R_ME_spice) <= 1e-12

    @test I_r_ME == frameconvert(R_r_BE, str2et(epoch_tdb), EM_BCR(), MJ2K())

    ###### Test Moon-centered J2000 to Rotating
    R_r_BE_rot = frameconvert(I_r_ME_spice, epoch_tdb, MJ2K(), EM_BCR())
    println("MJ2K -> Rotating: ", norm(R_r_BE_rot - R_r_BE))
    @test norm(R_r_BE_rot - R_r_BE) <= 1e12
    @test R_r_BE_rot == frameconvert(I_r_ME_spice, str2et(epoch_tdb), MJ2K(), EM_BCR())


    ###### Test Rotating to Earth-centered J2000
    I_r_EM = frameconvert(R_r_BM, epoch_tdb, EM_BCR(), EJ2K())
    I_R_EM = dimensionalize_state(AstroTOAST.em_cr3bp, I_r_EM)
    println("Rotating -> EJ2K: ", norm(I_R_EM - I_R_EM_spice))
    @test norm(I_R_EM - I_R_EM_spice) <= 1e-12
    @test I_r_EM == frameconvert(R_r_BM, str2et(epoch_tdb), EM_BCR(), EJ2K())

    ###### Test Earth-centered J2000 to Rotating
    R_r_BM_rot = frameconvert(I_r_EM_spice, epoch_tdb, EJ2K(), EM_BCR())
    println("EJ2K -> Rotating: ", norm(R_r_BM_rot - R_r_BM))
    @test norm(R_r_BM_rot - R_r_BM) <= 1e12
    @test R_r_BM_rot == frameconvert(I_r_EM_spice, str2et(epoch_tdb), EJ2K(), EM_BCR())

end
