# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                     CONSTRAINTS
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
"""
    Constraint{D}

Type Parameters
- `D`: dimension of the model
"""
abstract type Constraint{D} end

"""
    ContinuityConstraint

"""
struct ContinuityConstraint{D} <: Constraint{D}
    x1::FreeVariable
    x2::FreeVariable
    T::Float64
    model::DynamicalModel
end
