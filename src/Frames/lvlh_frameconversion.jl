using LinearAlgebra
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #
#                       FRAME CONVERSION BETWEEN LVLH AND EM_BCR
# -------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------- #

"""
`EM_BCR` <-> `EM_LVLH`
"""
function fc(target, chaser, f1::EM_BCR, f2::EM_LVLH)
    mu = mass_ratio(em_cr3bp)

    r_m_b_M = [1-mu,0,0] # Position of the moon wrt barycenter in the M frame

    # Target initial conditions (relative to the Moon)
    r0_M = target[1:3]-r_m_b_M # Initial r (target position wrt Moon)
    Mrdot0_M = target[4:6] # Initial rdot vector in the M frame
    target_dot = em_cr3bp(target) # Calculate rddot in the M frame
    Mrddot0_M = target_dot[4:6] # Initial rddot vector in the M frame

    # Define the LVLH frame (making sure to make target position wrt Moon)
    L_C_M = rot2lvlh(r0_M,Mrdot0_M)

    # h and r
    r0 = norm(r0_M) # Distance from Moon to target
    h0_M = cross(r0_M,Mrdot0_M) # Initial target angular momentum in the Moon frame
    h0 = norm(h0_M) # Initial target angular momentum scalar

    # Calculate wLM_L (from Franzini eq 31)
    wLMy = -h0/r0^2 # j component
    wLMz = -(r0/h0^2)*dot(h0_M,Mrddot0_M) # k component
    wLM_L = [0,wLMy,wLMz] # ang vel of L wrt M expressed in the L frame

    # Chaser initial conditions
    rho0_M = chaser[1:3]-target[1:3]
    rho0_L = L_C_M*rho0_M # Initial rho (r_chaser/target)
    rhodot0_M = chaser[4:6]-target[4:6] # rho dot = rcdot-rdot
    rhodot0_L = L_C_M*rhodot0_M - cross(wLM_L,rho0_L)  ## FIXED 02/15 to include cross product term

    return target, vcat(rho0_L,rhodot0_L)
end

function fc(target, chaser, f1::EM_LVLH, f2::EM_BCR)
    mu = mass_ratio(em_cr3bp)

    r_m_b_M = [1-mu,0,0] # Position of the moon wrt barycenter in the M frame

    # Unpack state
    r_M = target[1:3]-r_m_b_M # Target position wrt Moon -- r_t/m = r_t/o-r_m/o
    Mrdot_M = target[4:6] # Target velocity
    rho_L = chaser[1:3] # Eventually change to incorporate multiple chasers
    Lrhodot_L = chaser[4:6] # Chaser velocity [L frame]

    # Calculate acceleration of target -->cr3bp uses coordinates centered at barycenter
    Mtargdot_M = em_cr3bp(target) # S(1:6) is already wrt barycenter
    Mrddot_M = Mtargdot_M[4:6] # Target acceleration

    # Calculate DCM: L_C_M --> vec_L = L_C_M*vec_M
    L_C_M = rot2lvlh(r_M,Mrdot_M)

    # Calculate chaser position wrt barycenter in the M frame
    rCB_M = (L_C_M')*rho_L + r_M + r_m_b_M

    # Norms of states
    r = norm(r_M) # radius of target

    # Calculate angular momentum
    h_M = cross(r_M,Mrdot_M) # Angular momentum [M frame]
    h = norm(h_M) # Norm of angular momentum

    # Calculate wLM_L (from Franzini eq 31)
    wLMy = -h/r^2 # j component
    wLMz = -(r/h^2)*dot(h_M,Mrddot_M) # k component
    wLM_L = [0,wLMy,wLMz] # ang vel of L wrt M expressed in the L frame

    # Calculate chaser velocity wrt the M frame
    Mrcdot_M = (L_C_M')*(Lrhodot_L + cross(wLM_L,rho_L)) + Mrdot_M

    # Package output
    # S_M(i,:) = [S(i,1:6), rCB_M', Mrcdot_M']
    return target, vcat(rCB_M, Mrcdot_M)
end

"""
    rot2lvlh(r_M, Mrdot_M)

Calculate the rotation matrix necessary to transform from the Moon centered rotating frame into
the LVLH frame (position only). Inputs are position of the target with respect to the 

Returns the matrix `L_C_M`, where multiplying a vector
expressed in 3bp rotating frame coordinates with the 
matrix `L_C_M` yields a vector expressed in LVLH coords
The rows of the `L_C_M` matrix are the unit vectors of the L frame

INPUTS:
  `r_M` - vector position of the target relative to the Moon, expressed in M frame
  `Mrdot_M` - vector velocity of the target, taken wrt M, expressed in M

OUTPUTS:
  `L_C_M` - DCM to convert a 3 vector from moon frame to lvlh frame [3x3]
"""
function rot2lvlh(r_M, Mrdot_M)
    # Calculate angular momentum
    h_M = cross(r_M,Mrdot_M) # Angular momentum [M frame]
    h = norm(h_M) # Norm of angular momentum

    # Scalar position
    r = norm(r_M)

    # Calculate LVLH unit vectors
    jhat_M = -h_M/h # jhat as expressed in the moon frame
    khat_M = -r_M/r # khat as expressed in the moon frame
    ihat_M = cross(jhat_M,khat_M) # ihat as expressed in the moon frame

    # Calculate DCM: L_C_M --> vec_L = L_C_M*vec_M
    return L_C_M = [ihat_M';jhat_M';khat_M']
end
