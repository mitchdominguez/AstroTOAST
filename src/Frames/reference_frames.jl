# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                                  REFERENCE FRAME
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
    ReferenceFrame

Abstract type that will encompass different reference frames
"""
abstract type ReferenceFrame end

"""
    EM_BCR

Standard EM CR3BP rotating frame, centered at the Earth-Moon (EM) barycenter
"""
struct EM_BCR <: ReferenceFrame end

"""
    EM_ECR

EM CR3BP rotating frame, centered at the Earth
"""
struct EM_ECR <: ReferenceFrame end

"""
    EM_MCR

EM CR3BP rotating frame, centered at the Moon
"""
struct EM_MCR <: ReferenceFrame end

"""
    EM_ECI

Arbitrary Earth centered initial frame
"""
struct EM_ECI <: ReferenceFrame end

"""
    EM_MCI

Arbitrary moon centered initial frame
"""
struct EM_MCI <: ReferenceFrame end

"""
    EM_TCR

EM CR3BP rotating frame, centered on a target spacecraft
"""
struct EM_TCR <: ReferenceFrame end

"""
    EM_LVLH

LVLH frame defined wrt a target spacecraft in the EM CR3BP
"""
struct EM_LVLH <: ReferenceFrame end

"""
    EM_ICR

In-track, cross-track, radial frame defined wrt a target spacecraft in the EM CR3BP

The difference between EM_ICR and EM_LVLH is that inertial velocity is used
in the EM_ICR frame to calculate the angular momentum direction
"""
struct EM_ICR <: ReferenceFrame end



### Associated files
include("frame_conversion.jl")
