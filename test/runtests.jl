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
end
