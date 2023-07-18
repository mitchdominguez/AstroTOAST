using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test

# Calculate the 9:2 NRHO
#   This will define X1, X2, X3, T1, T2, T3, xv, among others
include("nrho92.jl") 

@testset "trajectory.jl" begin
    # Test that we actually have the 9:2 NRHO
    # @test Td == 1.511261560928471
    @test Td == 1.5091706485621377

    # Create some different sample trajectories
    traj_0 = Trajectory(model, tofullvector(X1), (0, tofullvector(T1)[1]))

    traj_1 = Trajectory(model, [X1], [T1])
    traj_2 = Trajectory(model, [X2], [T2])
    traj_3 = Trajectory(model, [X3], [T3])

    traj123 = Trajectory(model, [X1, X2, X3], [T1, T2, T3])

    # Test dimension
    @test dimension(traj_0) == dimension(traj_1) == dimension(traj_2) == 
            dimension(traj_3) == dimension(traj123) == 6

    # Test initial conditions
    @test x0(traj_0) == x0(traj_1) == x0(traj123)
    @test x0(traj_2) != x0(traj123)
    @test x0(traj123) == tofullvector(X1)

    # Test dynamical model
    @test dm(traj123) == Cr3bpModel(Bodies["Earth"],Bodies["Moon"])

    # Test tspan
    @test tspan(traj_0) == tspan(traj_1)
    @test tspan(traj_1) == [0, tofullvector(T1)[1]]
    @test tspan(traj_2) == [0, tofullvector(T2)[1]]
    @test tspan(traj_3) == [0, tofullvector(T3)[1]]
    @test tspan(traj123) ≈ [0, Td] atol=1e-12

    # Test overall time of flight
    @test tof(traj_1) == tofullvector(T1)[1]
    @test tof(traj_2) == tofullvector(T2)[1]
    @test tof(traj_3) == tofullvector(T3)[1]
    @test tof(traj123) ≈ Td atol=1e-12

    # Test solvec final results
    @test solvec(traj_0)[end][end] == solvec(traj_1)[end][end]
    @test solvec(traj_1)[end][end] ≈ x0(traj_2) atol=1e-12
    @test solvec(traj_2)[end][end] ≈ x0(traj_3) atol=1e-12
    @test solvec(traj123)[end][end] ≈ solvec(traj_3)[end][end] atol=1e-12

    # Test length
    @test length(traj_0) == length(traj_1) == length(traj_2) == length(traj_3) == 1
    @test length(traj123) == 3

    # Test append
    append!(traj_1, traj_2, traj_3)
    @test solvec(traj123)[end][end] ≈ solvec(traj_1)[end][end] atol=1e-12
    @test tspan(traj_1) ≈ [0, Td] atol=1e-12
    @test tof(traj_1) ≈ Td atol=1e-12
    @test length(traj_1) == 3 # Traj 1 gets appended to
    @test length(traj_2) == 1 # Traj 2 remains the same
    @test length(traj_3) == 1 # Traj 3 remains the same

    # Test append with compound trajectories
    ctraj = get_traj(nrho92)
    @test length(ctraj) == 3
    append!(ctraj, get_traj(nrho92))
    @test length(ctraj) == 6
    @test iscontinuous(ctraj) == true

    # Test getindex
    @test traj123[1].t[end] == traj_1[1].t[end]
    @test length(traj123[1:2]) == 2
    @test traj123[1:2][2].t[end] == traj123[2].t[end]

    # Test traj(time)
    @test traj123(0) == tofullvector(X1)
    @test traj123(tof(traj123)) ≈ tofullvector(X1) atol=1e-12
    @test traj123([0, tof(traj123)]) ≈ [tofullvector(X1), tofullvector(X1)] atol=1e-12

    # Test isperiodic
    @test isperiodic(traj123)
    @test isperiodic(traj123,1e-16) == false
    @test isperiodic(traj123,1e-5)


    
    
end
