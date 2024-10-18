using SparseArrays
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
    # evalmat = Array{Array{Float64},2}(undef,size(DFX(T)))
    # # evalmat = Array{AbstractArray,2}(undef,size(DFX(T)))

    # for j = 1:numels(X(T))
        # for i = 1:numels(FX(T))
            # # println(i)
            # # println(j)
            # # println(DFX(T)[i,j](FX(T)[i]))
            # evalmat[i,j] = DFX(T)[i,j](FX(T)[i])
        # end
    # end

    # # return evalmat

    # temp = Array{Array{Float64}}(undef,1,numels(X(T)))
    # for i = 1:size(evalmat,2)
        # # println(i)
        # # println(evalmat[1,i])
        # # println(evalmat[2,i])
        # # println(evalmat[3,i])
        # temp[1,i] = vcat(evalmat[:,i]...)
        # # println(size(temp[1,i]))
    # end 
    # outmat=hcat(temp...)[setdiff(1:end,removeinds(FX(T))), setdiff(1:end,removeinds(X(T)))]



    temp = zeros(full_length(FX(T)),full_length(X(T)))
    c_ind = 1

    for i = 1:numels(X(T))
        n = full_length(X(T)[i])
        r_ind = 1
        for j = 1:numels(FX(T))

            m = full_length(FX(T)[j])

            if !isa(DFX(T)[j,i], __NP)
                temp[r_ind:(r_ind+m-1),c_ind:(c_ind+n-1)] = DFX(T)[j,i](FX(T)[j])
            end

            r_ind += m
        end
        c_ind += n
    end

    outmat=temp[setdiff(1:end,removeinds(FX(T))), setdiff(1:end,removeinds(X(T)))]
end

"""
    evalSparseDFXMatrix

Evaluate a Sparse DF Matrix
"""
function construct_sparse_DF(T::Targeter)

    # temp = spzeros(Df,Dx)
    temp = spzeros(full_length(FX(T)),full_length(X(T)))
    # temp = zeros(Df,Dx)
    c_ind = 1

    for i = 1:numels(X(T))
        n = full_length(X(T)[i])
        r_ind = 1
        for j = 1:numels(FX(T))

            m = full_length(FX(T)[j])

            if !isa(DFX(T)[j,i], __NP)
                temp[r_ind:(r_ind+m-1),c_ind:(c_ind+n-1)] = DFX(T)[j,i](FX(T)[j])
            end

            r_ind += m
        end
        c_ind += n
    end

    outmat=temp[setdiff(1:end,removeinds(FX(T))), setdiff(1:end,removeinds(X(T)))]
end

function compute_invDF_times_F(fvec::AbstractVector, dfmat::AbstractMatrix; method=:fancy)
    if method == :backslash
        return dfmat\fvec
    elseif method == :fancy
        return dfmat'*((dfmat*dfmat')\fvec)
    else
        throw(ErrorException("Must specify inversion method!"))
    end
end

"""
    target(T::Targeter, maxiter::Int, tol::Float64, att::Float64)

Solve the targeting problem, specifying maximum number of iterations,
convergence tolerance, and attenuation factor
"""
function target(T::Targeter; maxiter=getmaxiter(T)::Int,
        tol=gettol(T)::Float64, att=1.0::Float64, debug=false, 
        eval_DF_func=evalDFXMatrix, inversion_method=:backslash,
        graceful_exit=false, infinity_norm_tol=nothing)

    # Save a copy of the initial value
    Xhist = Vector{XVector}()
    push!(Xhist, copy(X(T)))

    # Save the error history
    err = Vector{Float64}()
    push!(err, norm(FX(T)))
    if debug
        println("Iteration $(1): |F(X)| = $(err[1])")
        # println("\t X = $(tofullvector(X(T)))")
    end
    if err[1] < tol
        return (Xhist, err) 
    end

    i = 1
    cont = true
    while i<=maxiter && cont

        # Perform update
        # update!(X(T),tovector(X(T)) - att*(DFmat\tovector(FX(T))))
        update!(X(T),tovector(X(T)) - att*compute_invDF_times_F(tovector(FX(T)), eval_DF_func(T); method=inversion_method))

        # Update history
        push!(Xhist, copy(X(T)))

        # Evaluate stopping condition
        push!(err, norm(FX(T))) # update error history
        
        infty_err = Vector{Float64}()
        if err[end] < tol
            if isnothing(infinity_norm_tol)
                cont = false
            else
                push!(infty_err, infty_constraint_mags(FX(T)))
                if infty_err[end] < infinity_norm_tol
                    cont = false
                end
            end
        end
        i+=1

        if debug
            if isnothing(infinity_norm_tol) || isempty(infty_err)
                println("Iteration $(i): |F(X)| = $(err[i])")
            else
                println("Iteration $(i): |F(X)| = $(err[i])   -   |F(X)|∞ = $(infty_err[end])")
            end
        end
    end

    if i>maxiter
        if graceful_exit
            @warn "Exceeded maximum number of iterations ($(maxiter))! Returning Xhist and err..."
            return (Xhist, err)
        else
            error("maximum number of iterations reached")
        end
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
