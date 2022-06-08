# GENERATE THE 2:1 Halo
using AstroTOAST
using LinearAlgebra
using StaticArrays
using Roots
using Test

# Define the model
model = Cr3bpModel(Bodies["Earth"],Bodies["Moon"])

# Desired period
P_syn = 29.531*day2sec/dimensional_time(model) # Synodic period of the EM system wrt the Sun (ndim)
Td = (1/2)*P_syn # 9:2 NRHO

# Define some free variables and continuity constraints
X1 = FreeVariable("x1", [1.178980143388651, 0.0, -0.0429267707, 0.0, -0.16568789625925529, 0.0], includeinds=[1,3,5]) # Dhruv
T1 = FreeVariable("T1", 3.4003445182912633-1e-3, includeinds=[]) # Dhruv 
cc = ContinuityConstraint(X1,X1,T1,model,includeinds=[1,3,6])

#### THIS WORKS FOR FIXING Z, DON'T TOUCH IT ################################################################
# X1 = FreeVariable("x1", [1.178977343109523, 0, -0.042926770700000, 0, -0.165685370741522, 0], includeinds=[1,5]) # 2:1 Halo
# T1 = FreeVariable("T1", 3.39563395926481)
# cc = ContinuityConstraint(X1,X1,T1,model,includeinds=[1,3,6]) # constrain s.t we have a perpendicular crossing, must constrain y, xdot, zdot
# ################################################################################################################

# Define XVector
xv = XVector(X1, T1)

# Define FXVector
fx = FXVector(cc) # FX vector for rm in X, full FX

# Define Targeter
maxiter = 15
tol = 1e-12
targ = Targeter(xv, fx, maxiter, tol);

# Target initial PO
target(targ);
println(tofullvector(fx))

# Range of periods to target
tees = LinRange(tofullvector(T1)[1],Td,10)

Xhist = []
err = []

# Continue in position towards 2:1 Halo
for i = 1:length(tees)
    # Update tof
    update!(T1, [tees[i]])

    # Target new PO
    Xhist_temp, err_temp = target(targ)
    push!(Xhist,Xhist_temp)
    push!(err,err_temp)

    # println(tofullvector(fx))
end

halo21_traj = Trajectory(model, [X1], [T1])
halo21 = PeriodicOrbit(halo21_traj, "2:1 Halo", "L2 Halos")
