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
struct Targeter{Df,Dx}
    X::XVector
    FX::FXVector
    DFX::Matrix{Partial}
    maxiter::Int
    tol::Float64

    function Targeter(X::XVector, FX::FXVector, maxiter::Int, tol::Float64)

        # Generate DF Matrix Function array
        DFX = Matrix{Partial}(undef,numels(FX),numels(X))
        for j = 1:numels(X)
            for i = 1:numels(FX)
                DFX[i,j] = partials(FX[i],X[j])
            end
        end

        new{length(FX), length(X)}(X, FX, DFX, maxiter, tol)
    end
end

"""
    X(T::Targeter)

Return XVector
"""
X(T::Targeter) = T.X

"""
    FX(T::Targeter)

Return FXVector
"""
FX(T::Targeter) = T.FX

"""
    DFX(T::Targeter)

Return DFX Matrix
"""
DFX(T::Targeter) = T.DFX

"""
    getmaxiter(T::Targeter)

Return maxiter
"""
getmaxiter(T::Targeter) = T.maxiter

"""
    gettol(T::Targeter)

Return maxiter
"""
gettol(T::Targeter) = T.tol

"""
    evalDFXMatrix

Evaluate the DF Matrix
"""
function evalDFXMatrix(T::Targeter)
    evalmat = Array{Array{Float64}}(undef,size(DFX(T)))

    for j = 1:numels(X(T))
        for i = 1:numels(FX(T))
            # println(i)
            # println(j)
            # println(DFX(T)[i,j](FX(T)[i]))
            evalmat[i,j] = DFX(T)[i,j](FX(T)[i])
        end
    end

    temp = Array{Array{Float64}}(undef,1,numels(X(T)))
    for i = 1:size(evalmat,2)
        # println(i)
        # println(evalmat[1,i])
        # println(evalmat[2,i])
        # println(evalmat[3,i])
        temp[1,i] = vcat(evalmat[:,i]...)
        # println(size(temp[1,i]))
    end 
    outmat=hcat(temp...)[setdiff(1:end,removeinds(FX(T))), setdiff(1:end,removeinds(X(T)))]
end

"""
    target(T::Targeter, maxiter::Int, tol::Float64, att::Float64)

Solve the targeting problem, specifying maximum number of iterations,
convergence tolerance, and attenuation factor
"""
function target(T::Targeter; maxiter=getmaxiter(T)::Int, tol=gettol(T)::Float64, att=1.0::Float64)
    # Save a copy of the initial value
    Xhist = Vector{XVector}()
    push!(Xhist, copy(X(T)))

    # Save the error history
    err = Vector{Float64}()
    push!(err, norm(FX(T)))

    i = 1
    cont = true
    while i<=maxiter && cont
        # Evaluate DF Matrix
        DFmat = evalDFXMatrix(T)

        # Perform update
        update!(X(T),tovector(X(T)) - att*(DFmat\tovector(FX(T))))

        # Update history
        push!(Xhist, copy(X(T)))

        # Evaluate stopping condition
        push!(err, norm(FX(T))) # update error history
        # println(err[end]<tol)
        if err[end] < tol
            cont = false
        end
        i+=1
        # println(size(Xhist))
        # println(i)
    end

    if i>maxiter
        # println("Error history: $(err)")
        println("Error history")
        show(stdout, "text/plain", err)
        println("\n")
        error("maximum number of iterations reached")
    end

    return (Xhist, err)
end

"""
    Base.show

Overload the show operator to pretty print the Targeter to the console.
"""
function Base.show(io::IO, ::MIME"text/plain", targ::Targeter{Df,Dx}) where {Df,Dx}
    print(io, "Targeter:\n")
    print(io, "- Maximum iterations allowed: $(targ.maxiter)\n")
    print(io, "- Tolerance on ||F(X)||: $(targ.tol)\n")
    print(io, "\n=====\n")
    show(io, "text/plain", targ.X)
    print(io, "\n=====\n")
    show(io, "text/plain", targ.FX)
    print(io, "Summary: Targeter($(Df)x$(Dx))\n")
end
