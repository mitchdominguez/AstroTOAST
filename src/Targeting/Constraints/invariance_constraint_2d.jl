# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                  INVARIANCE CONSTRAINT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
    InvarianceConstraint2D

Invariance constraint for a quasiperiodic orbit on a 2-dimensional torus
"""
struct InvarianceConstraint2D{D} <: Constraint{D}
    u0::FreeVariable{D1,T} where {D1,T}
    xstar::FreeVariable{D2,T} where {D2,T} # TODO make this a vector?
    T::FreeVariable{D3,T} where {D3,T}
    ρ::FreeVariable{D4,T} where {D4,T}
    dm::DynamicalModel
    removeinds::Vector{Int}

    function InvarianceConstraint2D(u0, xstar, T, ρ, dm, rminds)
        # STUFF
        new()
    end
end

