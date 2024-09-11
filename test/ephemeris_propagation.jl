using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test

include("display.jl")

@testset "test_ephemeris_rotation.jl" begin

    ls_ricardo = 3.847479920112924e+05
    ts_ricardo = 3.756998590849912e+05


    x_embcr = [0.884916533745800,0,0,0,0.470613569321843,0]
    utcepoch = UTCEpoch(2025, 1, 1, 12, 0, 0.0)

    model = AstroTOAST.em_cr3bp

    hfem = HFEModel(AstroTOAST.Moon, dimensional_quantity_set(model), [AstroTOAST.Earth, AstroTOAST.Sun], utcepoch)

    I_R_ME_spice = ephemeris_state(;target=AstroTOAST.Earth, epoch=str2et(utcepoch), frame="J2000", observer=AstroTOAST.Moon)

    x_mj2k = frameconvert(x_embcr, utcepoch, EM_BCR(), MJ2K())

    x_mj2k_nadia = [-0.051220074685291,0.077314030343455,0.041924988067173,0.320729606562720,0.160572871492497,0.086751674764737]

    # println(x_embcr)
    # println(x_mj2k)
    # println(x_mj2k_nadia)
    # println(x_mj2k_nadia - x_mj2k)
    # println(I_R_ME_spice)
    
    acc_nadia = hfem(x_mj2k_nadia, 0.0)

    println("Acceleration: ", acc_nadia)

end
