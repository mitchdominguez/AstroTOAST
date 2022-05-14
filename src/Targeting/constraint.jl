# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                     CONSTRAINTS
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
using StaticArrays
"""
    Constraint{D}

Type Parameters
- `D`: dimension of the model
"""
abstract type Constraint{D} end

# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                       FX Vector
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

struct FXVector{D}
    Cs::Vector{Constraint} 

    # Constructor
    function FXVector(constraints...)
        # Initialize array
        # fv_array = Vector{FreeVariable}(undef,Base.length(constraints))
        fx_array = Vector{Constraint}()
        for i in eachindex(constraints)
            if ~(typeof(constraints[i])<:Constraint)
                throw(MethodError(FXVector, constraints[i]))
            end
            push!(fx_array, constraints[i])
        end
        # new{Base.length(constraints)}(fv_array)
        new{Base.length(fx_array)}(fx_array)
    end
end
