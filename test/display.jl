function ignore_zero(vee::AbstractVecOrMat, realtol = 1e-12, imagtol=1e-6)
    v = Vector(copy(vee))
    real_v = real(v)
    imag_v = imag(v)


    for i = 1:length(v)
        if abs(real_v[i])<=realtol
            real_v[i] = 0.
        end
        if abs(imag_v[i])<=imagtol
            imag_v[i] = 0.
        end
    end

    maxval, maxind = findmax(x->abs(x), imag(v))
    
    if maxval < imagtol
        return ignore_imag(complex.(real_v, imag_v))
    else
        return complex.(real_v, imag_v)
    end
end
