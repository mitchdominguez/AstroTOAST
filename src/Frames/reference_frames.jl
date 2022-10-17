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
name(::EM_BCR) = "EM Rotating"

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
name(::EM_TCR) = "Target-Centered Rotating"

"""
    EM_LVLH

LVLH frame defined wrt a target spacecraft in the EM CR3BP. 

This frame will always be reported with `target` in `EM_BCR` and `chaser` in
the LVLH frame
"""
struct EM_LVLH <: ReferenceFrame end
dm(::EM_LVLH) = em_cr3bp
name(::EM_LVLH) = "LVLH"

"""
    EM_VCR

In-track (̂iᵥ), cross-track (̂jᵥ), radial (̂kᵥ) frame defined wrt a target spacecraft in the EM CR3BP

The unit vectors defining the `EM_VCR` frame are defined as follows:

1) ̂iᵥ= velocity direction
2) ̂jᵥ= negative angular momentum direction (-̂h)
3) ̂kᵥ = completes the triad (îᵥ x ĵᵥ)

This frame will always be reported with `target` in `EM_BCR` and `chaser` in
the VCR frame
"""
struct EM_VCR <: ReferenceFrame end

"""
    EM_VNC

V(elocity), N(ormal), C(ross) frame, centered on a target in the EM CR3BP. The
normal is computed as the cross product of position wrt the Moon with the
rotating velocity
"""
struct EM_VNC <: ReferenceFrame end


# List of relative frames
relative_frames = [EM_TCR(), EM_LVLH(), EM_VCR(), EM_VNC()]

"""
    isrelativeframe(f::ReferenceFrame)

Return true if `f` is a relative frame
"""
isrelativeframe(f::ReferenceFrame) = f in relative_frames


### Associated files
include("frame_conversion.jl")
