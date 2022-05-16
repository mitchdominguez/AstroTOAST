# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                     Targeter
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
    Targeter{Dx, Df}

Targeter object, with type Dx denoting the length of the XVector, and Df denoting
the length of the FXVector vector
"""
struct Targeter{Dx,Df}
    X::XVector
    FX::FXVector
    maxiter::Int
    tol::Float64

    function Targeter(X::XVector, FX::FXVector, maxiter::Int, tol::Float64)
        new{length(X), length(FX)}(X, FX, maxiter, tol)
    end
end

"""
    Base.show

Overload the show operator to pretty print the Targeter to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", targ::Targeter)
    print(io, "Targeter:\n")
    print(io, "- Maximum iterations allowed: $(targ.maxiter)\n")
    print(io, "- Tolerance on ||F(X)||: $(targ.tol)\n")
    print(io, "\n=====\n")
    show(io, "text/plain", targ.X)
    print(io, "\n=====\n")
    show(io, "text/plain", targ.FX)
end
