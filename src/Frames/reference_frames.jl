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

dm(::EM_BCR) = em_cr3bp

"""
    EM_ECR

EM CR3BP rotating frame, centered at the Earth
"""
struct EM_ECR <: ReferenceFrame end
dm(::EM_ECR) = em_cr3bp

"""
    EM_MCR

EM CR3BP rotating frame, centered at the Moon
"""
struct EM_MCR <: ReferenceFrame end
dm(::EM_MCR) = em_cr3bp

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

EM CR3BP rotating frame, centered on a target spacecraft.

This frame will always be reported with `target` in `EM_BCR` and `chaser` in
the target centered frame that rotates with the EM CR3BP
"""
struct EM_TCR <: ReferenceFrame end
dm(::EM_TCR) = em_cr3bp

"""
    EM_LVLH

LVLH frame defined wrt a target spacecraft in the EM CR3BP. 

This frame will always be reported with `target` in `EM_BCR` and `chaser` in
the LVLH frame
"""
struct EM_LVLH <: ReferenceFrame end
dm(::EM_LVLH) = em_cr3bp

"""
    EM_ICR

In-track, cross-track, radial frame defined wrt a target spacecraft in the EM CR3BP

The difference between `EM_ICR` and `EM_LVLH` is that inertial velocity is used in
the `EM_ICR` frame to calculate the angular momentum direction.  

This frame will always be reported with `target` in `EM_BCR` and `chaser` in
the target centered ICR frame.
"""
struct EM_TCICR <: ReferenceFrame end



### Associated files
include("frame_conversion.jl")
