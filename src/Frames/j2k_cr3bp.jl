using SPICE
using LinearAlgebra
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                       FRAME CONVERSION BETWEEN EM_BCR, EJ2K, MJ2K
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
Construct DCM for EM_BCR -> J2000
"""
function construct_dcm_embcr_to_j2k(epoch::AbstractEpoch)

    #####
    # Lowercase r,v indicate nondimensional
    # Uppercase R,V indicate dimensional (km, s)
    #####

    ## Characteristic quantities for EM_BCR
    lstar = dimensional_length(em_cr3bp)
    tstar = dimensional_time(em_cr3bp)
    mu = em_cr3bp.μ
    GMstar = Moon.gravitational_parameter + Earth.gravitational_parameter

    ## State of Earth relative to the Moon
    x_EM_kms = ephemeris_state(;target=Moon, epoch=str2et(epoch), frame="J2000", observer=Earth)

    ## Instaneous characteristic quantities
    lstar_inst = norm(x_EM_kms[1:3]) # instantaneous lstar
    tstar_inst = sqrt((lstar_inst^3)/GMstar)

    ## Unit vectors of instantaneous EM rotating frame
    xhat_inst = normalize(x_EM_kms[1:3])
    zhat_inst = normalize(cross(x_EM_kms[1:3],x_EM_kms[4:6]))
    yhat_inst = normalize(cross(zhat_inst,xhat_inst))

    ## Instantaneous angular momentum vector
    H_inst = cross(x_EM_kms[1:3], x_EM_kms[4:6])

    ## DCM for converting position
    iCr = hcat(xhat_inst, yhat_inst, zhat_inst)

    ## Instantaneous radial velocity of Moon relative to Earth
    Vr_inst = dot(xhat_inst, x_EM_kms[4:6])

    ## Instantaneous dimensional anuglar velocity
    iWr = H_inst/(lstar_inst^2)

    ## Construct DCM rotating -> inertial
    DCM_11 = (lstar_inst/lstar)*iCr
    DCM_12 = zeros(3,3)
    DCM_21 = (tstar/lstar)*(Vr_inst*iCr + lstar_inst*(crs(iWr)*iCr))
    DCM_22 = (lstar_inst/lstar)*(tstar/tstar_inst)*iCr

    iCr_full = vcat(hcat(DCM_11, DCM_12), hcat(DCM_21, DCM_22))

    return iCr_full
end
function construct_dcm_embcr_to_j2k(epoch_et::Float64)

    #####
    # Lowercase r,v indicate nondimensional
    # Uppercase R,V indicate dimensional (km, s)
    #####

    ## Characteristic quantities for EM_BCR
    lstar = dimensional_length(em_cr3bp)
    tstar = dimensional_time(em_cr3bp)
    mu = em_cr3bp.μ
    GMstar = Moon.gravitational_parameter + Earth.gravitational_parameter

    ## State of Earth relative to the Moon
    x_EM_kms = ephemeris_state(;target=Moon, epoch=epoch_et, frame="J2000", observer=Earth)

    ## Instaneous characteristic quantities
    lstar_inst = norm(x_EM_kms[1:3]) # instantaneous lstar
    tstar_inst = sqrt((lstar_inst^3)/GMstar)

    ## Unit vectors of instantaneous EM rotating frame
    xhat_inst = normalize(x_EM_kms[1:3])
    zhat_inst = normalize(cross(x_EM_kms[1:3],x_EM_kms[4:6]))
    yhat_inst = normalize(cross(zhat_inst,xhat_inst))

    ## Instantaneous angular momentum vector
    H_inst = cross(x_EM_kms[1:3], x_EM_kms[4:6])

    ## DCM for converting position
    iCr = hcat(xhat_inst, yhat_inst, zhat_inst)

    ## Instantaneous radial velocity of Moon relative to Earth
    Vr_inst = dot(xhat_inst, x_EM_kms[4:6])

    ## Instantaneous dimensional anuglar velocity
    iWr = H_inst/(lstar_inst^2)

    ## Construct DCM rotating -> inertial
    DCM_11 = (lstar_inst/lstar)*iCr
    DCM_12 = zeros(3,3)
    DCM_21 = (tstar/lstar)*(Vr_inst*iCr + lstar_inst*(crs(iWr)*iCr))
    DCM_22 = (lstar_inst/lstar)*(tstar/tstar_inst)*iCr

    iCr_full = vcat(hcat(DCM_11, DCM_12), hcat(DCM_21, DCM_22))

    return iCr_full
end

"""
`EM_BCR` <-> `MJ2K`
"""
function fc(R_x_bs, epoch::Union{Float64,AbstractEpoch}, f1::EM_BCR, f2::MJ2K)
    
    iCr_full = construct_dcm_embcr_to_j2k(epoch)

    R_x_bm = [1-em_cr3bp.μ, 0, 0, 0, 0, 0]

    return iCr_full*(R_x_bs - R_x_bm)
end

function fc(I_x_ms, epoch::Union{Float64,AbstractEpoch}, f1::MJ2K, f2::EM_BCR)
    
    iCr_full = construct_dcm_embcr_to_j2k(epoch)

    R_x_bm = [1-em_cr3bp.μ, 0, 0, 0, 0, 0]

    return inv(iCr_full)*I_x_ms + R_x_bm
end



"""
`EM_BCR` <-> `EJ2K`
"""
function fc(R_x_bs, epoch::Union{Float64,AbstractEpoch}, f1::EM_BCR, f2::EJ2K)

    iCr_full = construct_dcm_embcr_to_j2k(epoch)

    R_x_be = [-em_cr3bp.μ, 0, 0, 0, 0, 0]

    return iCr_full*(R_x_bs - R_x_be)
end

function fc(I_x_es, epoch::Union{Float64,AbstractEpoch}, f1::EJ2K, f2::EM_BCR)
    
    iCr_full = construct_dcm_embcr_to_j2k(epoch)

    R_x_be = [-em_cr3bp.μ, 0, 0, 0, 0, 0]

    return inv(iCr_full)*I_x_es + R_x_be
end


"""
`MJ2K` <-> `EJ2K`
"""
function fc(x_MS, epoch::AbstractEpoch, f1::MJ2K, f2::EJ2K)
    X_ME = ephemeris_state(;target=Earth, epoch=str2et(epoch), frame="J2000", observer=Moon)
    x_ME = nondimensionalize_state(em_cr3bp, X_ME)

    return x_MS - x_ME
end

function fc(x_ES, epoch::AbstractEpoch, f1::EJ2K, f2::MJ2K)
    X_ME = ephemeris_state(;target=Earth, epoch=str2et(epoch), frame="J2000", observer=Moon)
    x_ME = nondimensionalize_state(em_cr3bp, X_ME)

    return x_ES + x_ME
end


