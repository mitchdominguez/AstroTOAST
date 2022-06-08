# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                JACOBI CONSTANT CONSTRAINT
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
    JacobiConstraint

Applies the constraint `JCavg` - `JCd` = 0, where `JCavg` is the average of Jacobi
constant values for the FreeVariables contained within this constraint, and `JCd`
is the desired Jacobi constant.
"""
struct JacobiConstraint{D} <: Constraint{D}
    X::FreeVariable
    JCd::Real
    model::DynamicalModel
    removeinds::Vector{Int}

    # Constructor
    function JacobiConstraint(model::Cr3bpModel, JCd, X::FreeVariable)
        μ = mass_ratio(model)
        dim = dimension(model)

        # Check that Td is nonnegative
        if JCd < 0
            throw(InvalidStateException("JCd must be nonnegative", :JCd))
        end


        if full_length(X)%dim != 0
            throw(DimensionMismatch("All FreeVariables passed into the
                                    TOFConstraint constructor must be
                                    of full_length = a multiple of $(dim)"))
        end # Check X length


        new{1}(X, JCd, model, Vector{Int}([]))
    end

    JacobiConstraint() = new{1}()
end

"""
    jcd(jc::JacobiConstraint)

Return the desired Jacobi constant
"""
jcd(jc::JacobiConstraint) = jc.JCd

"""
    xvar(jc::JacobiConstraint)

Return FreeVariable
"""
xvar(jc::JacobiConstraint) = jc.X

"""
    dm(jc::JacobiConstraint)

Return the dynamical model of the Jacobi constant constraint
"""
dm(jc::JacobiConstraint) = jc.model

"""
    full_length(jc::JacobiConstraint)

Return the full length of the JacobiConstraint
"""
full_length(jc::JacobiConstraint) = length(jc)

"""
    evalconstraint(jc::JacobiConstraint)

Evaluate the Jacobi constraint
"""
function evalconstraint(jc::JacobiConstraint)
    JCd = jcd(jc)
    X = xvar(jc)
    model = dm(jc)
    D = dimension(model)

    N = Int(full_length(X)/D)
    jcvec = Vector{Float64}(undef, N)
    for i = 1:N
        jcvec[i] = jacobi_constant(model, X[D*i-(D-1):D*i])
    end
    jcv = vcat(jcvec...)

    return sum(jcv)/length(jcv) - JCd
end

"""
    partials(jc::JacobiConstraint, fv::FreeVariable)

Return the matrix of partial derivatives for the partial of the constraint with
respect to the given free variable
"""
function partials(jc::JacobiConstraint, fv::FreeVariable{D,T}) where {D,T}
    if fv==xvar(jc)
        return __dJC_dX{full_length(fv)}()
    else
        # No partial
        # return __NP{D}()
        return __NP{full_length(fv)}()
    end
end

"""
    __dJC_dX

Partial of the time of flight constraint with respect to Ti, which
is a FreeVariable included in the TOF constraint
"""
struct __dJC_dX{D} <: Partial{D} end
function (::__dJC_dX{C})(jc::JacobiConstraint{R}) where {R,C}
    X = xvar(jc)
    model = dm(jc)
    D = dimension(model)
    N = Int(full_length(X)/D)

    outmat = zeros(R,C)

    for i = 1:N
        r = D*i-(D-1):D*i
        outmat[r] = 2*vcat(pseudopotential_gradient(model, X[r]), -X[r[end-2:end]])
    end

    return outmat
end