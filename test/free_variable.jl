using AstroTOAST
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

    # length
    @test length(X0) == 1
    @test length(X1) == 3
    @test length(X2) == 6
    @test length(X3) == 9

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


    # XVector construction
    @test_throws MethodError XVector{3}([X1, X2, X3])
    @test_throws MethodError xv2 = XVector(X1, 5)
    xv = XVector(X1, X2, X3)

    # numels
    @test numels(xv) == 3

    # length
    @test length(xv) == 18

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

end
