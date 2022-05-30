using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test

# Calculate the 9:2 NRHO
#   This will define nrho_traj
include("nrho92.jl") 

@testset "trajectory.jl" begin
    # Test that we actually have the 9:2 NRHO
    @test tof(nrho92_traj) == 1.509170648562138
    
    # Some sample trajectories
    traj1 = Trajectory(model, [X1], [T1])
    
    # Constructor
    nrho92 = PeriodicOrbit(nrho92_traj, "9:2 NRHO", "L2 Halos")
    @test_throws ErrorException PeriodicOrbit(traj1)

    # traj
    @test tof(traj(nrho92)) == tof(nrho92_traj)
    @test x0(traj(nrho92)) == x0(nrho92_traj)
    @test traj(nrho92) != nrho92_traj

    # dimension
    @test dimension(nrho92) == 6

    # x0
    @test x0(nrho92) == x0(nrho92_traj)

    # dm
    @test dm(nrho92) == dm(nrho92_traj) == model

    # period
    @test period(nrho92) == tof(nrho92_traj)

    # JC
    @test jacobi_constant(nrho92) == 3.04664708021343

    # Test getting periodic orbit state at a specific time
    @test nrho92(0) == nrho92_traj(0) == tofullvector(X1)
    @test nrho92(1) == nrho92_traj(1)
    @test nrho92(period(nrho92)) ≈ tofullvector(X1) atol=1e-12
    @test nrho92(period(nrho92)) ≈ nrho92_traj(tof(nrho92_traj)) atol=1e-12

    # name, family
    @test name(nrho92) == "9:2 NRHO"
    @test family(nrho92) == "L2 Halos"

    # Monodromy matrix
    M = [-0.894376     -0.385543       8.64367     0.000996841     0.0853912   -0.00366975;
         520.949      -112.527      -7064.82       -0.0853912     -68.6522      -1.10901;
         -13.051         5.11751      113.958      -0.00366975      1.10901      0.0494274;
         5412.58      -1556.19      -67611.3        -1.06516      -658.253      -15.269;
         1332.93       -529.117     -11581.0         0.383549     -112.698       -5.11017;
         -53481.7       11563.7      7.27032e5      8.64367       7064.82       113.958]
    @test monodromy(nrho92)[5] ≈ M[5] atol=1e-2 # Printing M means losing lots of precision

    # Eigenvalues
    @test eigvals(nrho92) != eigvals(monodromy(nrho92)) # Different order
    @test isempty(setdiff(eigvals(monodromy(nrho92)), eigvals(nrho92))) # All same elements

    # Eigenvectors
    @test size(eigvecs(nrho92)) == (6,)
    @test size(eigvecs(monodromy(nrho92))) == (6,6)

    # Test that eigenvalues, eigenvectors are paired properly
    for i = 1:6
        @test norm(monodromy(nrho92)*eigvecs(nrho92)[i] - eigvecs(nrho92)[i]*eigvals(nrho92)[i]) ≈ 0 atol=1e-12
    end

    # classify_eigs
    @test classify_eigs(nrho92)[1] == [1] # unstable
    @test classify_eigs(nrho92)[2] == [5,6] # center 
    @test classify_eigs(nrho92)[3] == [3,4] # unit
    @test classify_eigs(nrho92)[4] == [2] # stable

    # Obtaining different classes of eigenvalues, vectors
    u_eig = unstable_eigs(nrho92)
    @test u_eig[1][1] == eigvals(nrho92)[1]
    @test u_eig[2][1] == eigvecs(nrho92)[1]

    c_eig = center_eigs(nrho92)
    @test c_eig[1][1] == eigvals(nrho92)[5]
    @test c_eig[1][2] == eigvals(nrho92)[6]
    @test c_eig[2][1] == eigvecs(nrho92)[5]
    @test c_eig[2][2] == eigvecs(nrho92)[6]
    
    unit_eig = unit_eigs(nrho92)
    @test unit_eig[1][1] == eigvals(nrho92)[3]
    @test unit_eig[1][2] == eigvals(nrho92)[4]
    @test unit_eig[2][1] == eigvecs(nrho92)[3]
    @test unit_eig[2][2] == eigvecs(nrho92)[4]

    s_eig = stable_eigs(nrho92)
    @test s_eig[1][1] == eigvals(nrho92)[2]
    @test s_eig[2][1] == eigvecs(nrho92)[2]

    # Stability index
    @test stability_index(nrho92)[2] == 1.3187408624058536
    @test stability_index(nrho92)[1] == [-1.3187408624058536, 1.0000000012691992, 0.6844867717750143]
    
end
