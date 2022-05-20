using AstroTOAST
using StaticArrays
using Test


@testset "free_variable.jl" begin
    # constructor
    @test_throws MethodError FreeVariable{1,Int64}("x",[1])
    @test_throws MethodError FreeVariable(4,[1])

    x1vec = [1,2,3]
    x2vec = [1.0,2.0,3.0,4.0,5.0,6.0]
    x3vec = [1,2,3,4,5,6,7,8,9]

    # Construct some free variables
    X0 = FreeVariable("x0",1)
    X1 = FreeVariable("x1",x1vec)
    X2 = FreeVariable("x2",[1.0,2.0,3.0,4.0,5.0,6.0])
    X3 = FreeVariable("x3",[1,2,3,4,5,6,7,8,9])
    X4 = FreeVariable("x4",[0.72, 0, 0.71, 0, 0.18, 0],[2,4,6])
    X5 = FreeVariable("x5",[0.72, -1, 0.71, 0.39, 0.18, 0.46],3)
    T = FreeVariable("T",2.0,1)

    # removeinds
    @test removeinds(X0) == []
    @test removeinds(X4) == [2,4,6]

    # active
    @test active(X1) == true
    @test active(T) == false
    @test active(X4) == true

    # length
    @test length(X0) == 1
    @test length(X1) == 3
    @test length(X2) == 6
    @test length(X3) == 9
    @test length(X4) == 3
    @test full_length(X4) == 6
    # @test length(X4, true) == 3

    # element type
    @test eltype(X1) == Int64
    @test eltype(X2) == Float64
    @test eltype(X3) == Int64

    # name
    @test name(X1) == "x1"
    @test name(X2) == "x2"
    @test name(X3) == "x3"

    # value
    @test value(X1) == [1,2,3]
    @test value(X2) == [1.0,2.0,3.0,4.0,5.0,6.0]
    @test value(X3) == [1,2,3,4,5,6,7,8,9]

    # tovector
    @test tovector(X1) == [1,2,3]
    @test tovector(X2) == [1.0,2.0,3.0,4.0,5.0,6.0]
    @test tovector(X3) == [1,2,3,4,5,6,7,8,9]
    @test tovector(X4) == [0.72, 0.71, 0.18]
    @test tofullvector(X4) == [0.72, 0, 0.71, 0, 0.18, 0]
    @test tofullsvector(X4) == @SVector [0.72, 0, 0.71, 0, 0.18, 0]

    # test that changing the value of x1vec does not change X1
    @test value(X1) == x1vec
    x1vec[2] = 5000
    @test value(X1) != x1vec

    # iterate 
    @test iterate(X1) == (1,2)
    @test iterate(X2,4) == (4.0,5)
    @test iterate(X3,12) == nothing

    # getindex
    @test X1[2] == 2
    @test getindex(X2,5) == 5.0
    @test_throws BoundsError X2[10]
    @test X3[8] == 8
    @test X3[2:5] == [2,3,4,5]
    @test_throws BoundsError X2[5:10]

    fullvec = [0.72, 0, 0.71, 0, 0.18, 0]
    for i = 1:length(fullvec)
        @test X4[i] == fullvec[i]
    end


    # XVector construction
    @test_throws MethodError XVector{3}([X1, X2, X3])
    @test_throws MethodError xv2 = XVector(X1, 5)
    xv = XVector(X1, X2, X3)
    xvf = XVector(X1, X2, X3, T)
    xvrm = XVector(X1, X4) # test removing free variables
    xvrm2 = XVector(X4, X5) # test removing free variables

    # removeinds
    @test removeinds(xvrm2) == [2,4,6,9]

    # numels
    @test numels(xv) == 3
    @test numels(xvf) == 3 # test that inactive free variables are not added
    @test numels(xvrm) == 2

    # length
    @test length(xv) == 18
    @test length(xvf) == 18
    @test length(xvrm) == 6 # The removed free variables are counted
    @test length(xvrm2) == 8 # The removed free variables are counted
    @test full_length(xvrm2) == 12
    # @test length(xvrm,true) == 6 # The removed free variables should not be counted

    # iterate
    @test iterate(xv,1) == (X1, 2)
    @test iterate(xv,2) == (X2, 3)
    @test iterate(xv,3) == (X3, 4)
    @test iterate(xv,4) == nothing

    # getindex
    @test getindex(xv,1) == X1 
    @test xv[1] == X1
    @test xv[1] != FreeVariable("x1",[1;2;3])
    @test xv[1:2] != [X1, X2]

    # tovector
    @test tovector(xv) == [1;2;3;x2vec;x3vec]
    # @test tovector(xv, true) == [1;2;3;x2vec;x3vec] # no elements to remove
    @test tofullvector(xvrm) == [1.0, 2.0, 3.0 ,0.72, 0, 0.71, 0, 0.18, 0]
    @test tovector(xvrm) == [1.0, 2.0, 3.0 ,0.72, 0.71, 0.18] # remove elements
    @test tofullsvector(xvrm) == @SVector [1.0, 2.0, 3.0 ,0.72, 0, 0.71, 0, 0.18, 0]
    @test tosvector(xvrm) == @SVector [1.0, 2.0, 3.0 ,0.72, 0.71, 0.18] # remove elements


    ######################################################

    # setindex! and update!
    @test value(X1) == [1,2,3]
    X1[2] = 17
    @test value(X1) == [1,17,3]
    @test_throws BoundsError X1[5]=99
    @test_throws MethodError X1[2]=4.5
    update!(X1, [9,99,999])
    @test value(X1) == [9, 99, 999]
    @test value(xv[1]) == [9, 99, 999]

    # update
    X1_update = update(X1, [22, 55, 77])
    @test value(X1) == [9, 99, 999]

    # copy
    X1_copy = copy(X1)
    @test X1_copy != X1
    @test name(X1_copy) == name(X1)
    @test value(X1_copy) == value(X1)
    
    # copy (XVector)
    xvc = copy(xv)
    @test xvc ≢ xv
    @test tovector(xvc) == tovector(xv)

    # setindex! and update! (XVector)
    X1f = FreeVariable("x1",[1.0, 2.0, 3.0])
    X2f = FreeVariable("x2",[1.3, 2.3, 3.3])
    X3f = FreeVariable("x3",[1.6, 2.6, 3.6])
    xvf = XVector(X1f, X2f, X3f) # all floats version
    xvf2 = XVector(copy(X3f), copy(X2f), copy(X1f))
    @test_throws DimensionMismatch xvc[2] = X2f
    update!(xvf, tovector(xvf).+1)
    @test value(xvf[1]) == [2.0, 3.0, 4.0]
    xvfc = update(xvf, tovector(xvf).+1)
    @test value(xvf[1]) == [2.0, 3.0, 4.0]
    @test value(xvfc[1]) == [3.0, 4.0, 5.0]
    @test_throws DimensionMismatch  xvfc = update(xvf, convert(Vector{Float64},[4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,0,0,0]))
    @test_throws BoundsError xvf[5] = X1f

    xvf_copy = copy(xvf)
    @test xvf_copy != xvf != xvf2
    @test xvf_copy[1] != xvf[1]
    @test xvf_copy[2] != xvf[2]
    @test xvf_copy[3] != xvf[3]
    @test tovector(xvf[1]) == [2.0, 3, 4]
    @test tovector(xvf_copy[1]) == [2.0, 3, 4]
    @test tovector(xvf2[1]) == [1.6, 2.6, 3.6]

    update!(xvf, xvf2)
    @test xvf_copy != xvf != xvf2
    @test tovector(xvf[1]) == [1.6, 2.6, 3.6]
    @test tovector(xvf_copy[1]) == [2.0, 3, 4]
    @test tovector(xvf2[1]) == [1.6, 2.6, 3.6]

    update!(xvf, xvf_copy)
    @test xvf_copy != xvf != xvf2
    @test tovector(xvf[1]) == [2.0, 3, 4]
    @test tovector(xvf_copy[1]) == [2.0, 3, 4]
    @test tovector(xvf2[1]) == [1.6, 2.6, 3.6]


    # Test update! with different length update vectors
    newvec = [1,2,3,4,5,6,7,8,9,10,11,12.0] # update full vector

    update!(xvrm2, newvec)
    @test deleteat!(copy(newvec), removeinds(xvrm2)) == [1,3,5,7,8,10,11,12.0]
    @test tovector(xvrm2) == deleteat!(copy(newvec), removeinds(xvrm2))
    @test tofullvector(xvrm2) == copy(newvec)

    newvec_rm = tovector(xvrm2).+0.5 # exclude selected elements 
    update!(xvrm2, newvec_rm) # Pass in vector of length length(xvrm2) (as opposed to full_length)
    @test tofullvector(xvrm2) == [1.5,2,3.5,4,5.5,6,7.5,8.5,9,10.5,11.5,12.5]


end
