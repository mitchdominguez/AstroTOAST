using AstroTOAST
using Test

@testset "AstroTOAST.jl" begin
    # Write your tests here.
    include("free_variable.jl")
    include("continuity_constraint.jl")
    include("constraint.jl")
    include("targeter.jl")
    include("poms.jl")
    include("poss.jl")
    include("tofconstraint.jl")
    include("trajectory.jl")
    include("periodic_orbit.jl")
    include("qposs.jl")
    include("perpendicular_crossing.jl")
    include("lvlh_targeting.jl")

    load_only_default_kernels()
    include("ephemeris_rotation.jl")
    include("ephemeris_propagation.jl")
    include("ephemeris_continuity.jl")
end
