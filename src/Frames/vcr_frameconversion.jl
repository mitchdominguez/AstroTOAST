using LinearAlgebra
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                       FRAME CONVERSION BETWEEN EM_VCR AND EM_BCR
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
`EM_BCR` <-> `EM_VCR`
"""
function fc(target, chaser, f1::EM_BCR, f2::EM_VCR)

    # Initial conditions (relative to the Moon)
    x_M = frameconvert(target, EM_BCR(), EM_MCR())
    r_M = x_M[1:3]
    Mrdot_M = x_M[4:6]
    target_dot = em_cr3bp(target) # Calculate rddot in the M frame
    Mrddot_M = target_dot[4:6] # Initial rddot vector in the M frame

    # Define the VCR frame (making sure to make target position wrt Moon)
    V_C_M = rot2vcr(r_M,Mrdot_M)

    ihat_V = V_C_M[1,:]
    jhat_V = V_C_M[2,:]
    khat_V = V_C_M[3,:]

    # h and r
    v = norm(Mrdot_M) # Velocity magnitude
    h_M = cross(r_M,Mrdot_M) # Target angular momentum in the Moon frame
    Mhdot_M = cross(r_M, Mrddot_M) # Target angular momentum time derivative in the Moon frame
    h = norm(h_M) # Target angular momentum scalar

    # Calculate wLM_L (from Franzini eq 31)
    wVMx = -(1/h)*dot(Mhdot_M, khat_V) # ihat_V component
    wVMy = -(1/v)*dot(Mrddot_M, khat_V) # jhat_V component
    wVMz = (1/2)*( (1/v)*dot(Mrddot_M,jhat_V) + (1/h)*dot(Mhdot_M, ihat_V))
    wVM_V = [wVMx,wVMy,wVMz] # ang vel of V wrt M expressed in the V frame

    # Chaser initial conditions
    rho_M = chaser[1:3]-target[1:3]
    rho_V = V_C_M*rho_M # Initial rho (r_chaser/target)
    rhodot_M = chaser[4:6]-target[4:6] # rho dot = rcdot-rdot
    rhodot_V = V_C_M*rhodot_M - cross(wVM_V,rho_V)  ## FIXED 02/15/2021 to include cross product term

    return target, vcat(rho_V,rhodot_V)
end

function fc(target, chaser, f1::EM_VCR, f2::EM_BCR)

    # Initial conditions (relative to the Moon)
    x_M = frameconvert(target, EM_BCR(), EM_MCR())
    r_M = x_M[1:3]
    Mrdot_M = x_M[4:6]
    target_dot = em_cr3bp(target) # Calculate rddot in the M frame
    Mrddot_M = target_dot[4:6] # Initial rddot vector in the M frame
    rho_V = chaser[1:3] # Eventually change to incorporate multiple chasers
    Vrhodot_V = chaser[4:6] # Chaser velocity [L frame]

    # Define the VCR frame (making sure to make target position wrt Moon)
    V_C_M = rot2vcr(r_M,Mrdot_M)
    M_C_V = V_C_M'

    ihat_V = V_C_M[1,:]
    jhat_V = V_C_M[2,:]
    khat_V = V_C_M[3,:]

    # Calculate chaser position wrt barycenter in the M frame
    rCB_M = M_C_V*rho_V + target[1:3]

    # h and r
    v = norm(Mrdot_M) # Velocity magnitude
    h_M = cross(r_M,Mrdot_M) # Target angular momentum in the Moon frame
    Mhdot_M = cross(r_M, Mrddot_M) # Target angular momentum time derivative in the Moon frame
    h = norm(h_M) # Target angular momentum scalar

    # Calculate wLM_L (from Franzini eq 31)
    wVMx = -(1/h)*dot(Mhdot_M, khat_V) # ihat_V component
    wVMy = -(1/v)*dot(Mrddot_M, khat_V) # jhat_V component
    wVMz = (1/2)*( (1/v)*dot(Mrddot_M,jhat_V) + (1/h)*dot(Mhdot_M, ihat_V))
    wVM_V = [wVMx,wVMy,wVMz] # ang vel of V wrt M expressed in the V frame

    # Calculate chaser velocity wrt the M frame
    Mrcdot_M = M_C_V*(Vrhodot_V + cross(wVM_V,rho_V)) + Mrdot_M

    # Package output
    return target, vcat(rCB_M, Mrcdot_M)
end

"""
    rot2lvlh(r_M, Mrdot_M)

Calculate the rotation matrix necessary to transform from the Moon centered rotating frame into
the VCR frame (position only). Inputs are position of the target with respect to the 

Returns the matrix `V_C_M`, where multiplying a vector
expressed in 3bp rotating frame coordinates with the 
matrix `V_C_M` yields a vector expressed in VCR coords
The rows of the `V_C_M` matrix are the unit vectors of the VCR (V) frame

INPUTS:
  `r_M` - vector position of the target relative to the Moon, expressed in M frame
  `Mrdot_M` - vector velocity of the target, taken wrt M, expressed in M

OUTPUTS:
  `V_C_M` - DCM to convert a 3 vector from moon frame to lvlh frame [3x3]
"""
function rot2vcr(r_M, Mrdot_M)
    # Calculate angular momentum
    h_M = cross(r_M,Mrdot_M) # Angular momentum [M frame]
    h = norm(h_M) # Norm of angular momentum

    # Calculate VCR unit vectors
    ihat_V = Mrdot_M/norm(Mrdot_M)
    jhat_V = -h_M/h
    khat_V = cross(ihat_V, jhat_V)/norm(cross(ihat_V, jhat_V))

    # Calculate DCM: V_C_M --> vec_V = V_C_M*vec_M
    return V_C_M = [ihat_V';jhat_V';khat_V']
end
