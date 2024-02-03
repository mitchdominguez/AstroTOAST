using LinearAlgebra
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                       FRAME CONVERSION BETWEEN LVLH AND EM_BCR
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
`EM_MCR` <-> `EM_MCAI`
"""
# function fc(xtr, θ, f1::EM_MCR, f2::EM_MCAI)
function fc(xtr, θ, f1::R, f2::AI) where {R<:Union{EM_MCR, EM_ECR}, AI<:Union{EM_MCAI, EM_ECAI}}

    iCr = [[ cos(θ), -sin(θ), 0,      0,       0, 0];;
           [ sin(θ),  cos(θ), 0,      0,       0, 0];;
           [      0,       0, 1,      0,       0, 0];;
           [-sin(θ), -cos(θ), 0, cos(θ), -sin(θ), 0];;
           [ cos(θ), -sin(θ), 0, sin(θ),  cos(θ), 0];;
           [      0,       0, 0,      0,       0, 1]]'

    return iCr*xtr
end

# function fc(xti, θ, f1::EM_MCAI, f2::EM_MCR)
function fc(xti, θ, f1::AI, f2::R) where {R<:Union{EM_MCR, EM_ECR}, AI<:Union{EM_MCAI, EM_ECAI}}

    rCi = [[ cos(θ),  sin(θ), 0,       0,      0, 0];;
           [-sin(θ),  cos(θ), 0,       0,      0, 0];;
           [      0,       0, 1,       0,      0, 0];;
           [-sin(θ),  cos(θ), 0,  cos(θ), sin(θ), 0];;
           [-cos(θ), -sin(θ), 0, -sin(θ), cos(θ), 0];;
           [      0,       0, 0,       0,      0, 1]]'


    return rCi*xti
end

"""
`EM_MCAI` <-> `EM_RMCAI`
"""
function fc(target, chaser, epoch, f1::EM_MCAI, f2::EM_RMCAI)
    return (target, chaser-target)
end
function fc(target, chaser, epoch, f1::EM_RMCAI, f2::EM_MCAI)
    return (target, chaser+target)
end

"""
`EM_ECAI` <-> `EM_RECAI`
"""
function fc(target, chaser, epoch, f1::EM_ECAI, f2::EM_RECAI)
    return (target, chaser-target)
end
function fc(target, chaser, epoch, f1::EM_RECAI, f2::EM_ECAI)
    return (target, chaser+target)
end
