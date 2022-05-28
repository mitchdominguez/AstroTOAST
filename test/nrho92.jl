# GENERATE THE 9:2 NRHO
using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test

# Define the model
model = Cr3bpModel(Bodies["Earth"],Bodies["Moon"])

# Define some free variables
# X1 = FreeVariable("x1", [0.987, 0, 0.008, 0, 1.667, 0], [2,4,6]) # fix x0, y0, xdot0, zdot0
X1 = FreeVariable("x1", [0.987, 0, 0.008, 0, 1.667, 0]) # Don't fix any variables
X2 = FreeVariable("x2", [0.988, 0.023, -0.009,0.079, 0.488, -0.802])
X3 = FreeVariable("x3", [1.021, 0.011, -0.178,0.015, -0.100, -0.058])


T1 = FreeVariable("T1", 0.0239)
T2 = FreeVariable("T2", 0.6156)
T3 = FreeVariable("T3", 0.8718)

# Define XVector
xv = XVector(X1, X2, X3, T1, T2, T3)

# Define constraints

#   keep full state in X, remove some constraints
cc = Vector{ContinuityConstraint}()
push!(cc, ContinuityConstraint(X1, X2, T1, model))
push!(cc, ContinuityConstraint(X2, X3, T2, model))
push!(cc, ContinuityConstraint(X3, X1, T3, model, 5)) # Remove ydot constraint
# push!(cc, ContinuityConstraint(X3, X1, T3, model, [1,5])) # Remove ydot constraint

# Define FXVector
fx = FXVector(cc...) # FX vector for full X, rm in FX

# Define Targeter
maxiter = 20
tol = 1e-12
targ = Targeter(xv, fx, maxiter, tol);

# Target fx rc version
Xhist, err = target(targ);

println("After first targeting")
println(norm(tofullvector(fx)))
println(tofullvector(T1)[1] + tofullvector(T2)[1] + tofullvector(T3)[1])
println("---")

# Now that we have a periodic orbit, we need to reformulate the targeting problem
# to get a periodic orbit with ICs on the x-z plane and with a period
# that corresponds to the 9:2 NRHO

# Find the xz plane crossing and set the FreeVariables accordingly
po1 = Trajectory(model, [X1, X2, X3], [T1, T2, T3])
peri_time = find_zero(t -> po1(t)[2], (tof(po1)-0.2, tof(po1)), Roots.Bisection())
dt = peri_time-tof(po1)
update!(X1, po1(peri_time).*[1,0,1,0,1,0])
update!(X2, po1(tofullvector(T1)[1] + dt))
update!(X3, po1(tofullvector(T1)[1] + tofullvector(T2)[1] + dt))

println("After updating free variables")
println(norm(tofullvector(fx)))

X1_fix = FreeVariable("x1f", tofullvector(X1), [2,4,6])

# New XVector
xv2 = XVector(X1_fix, X2, X3, T1, T2, T3)

# Define time of flight constraint
Td = 1.511261560928471 # 9:2 NRHO
tofc = TOFConstraint(Td, T1, T2, T3)

# Define a new set of continuity constraints
cc2 = Vector{ContinuityConstraint}()
push!(cc2, ContinuityConstraint(X1_fix, X2, T1, model,[1,3,5]))
push!(cc2, ContinuityConstraint(X2, X3, T2, model))
push!(cc2, ContinuityConstraint(X3, X1_fix, T3, model)) # Remove ydot constraint
# push!(cc2, ContinuityConstraint(X3, X1, T3, model, [1,5])) # Remove ydot constraint

# Make a new FXVector to add in the TOF constraint
fx2 = FXVector(cc2..., tofc)

println("New constraint")
println(norm(tofullvector(fx2)))

# Target with tof constraint
targ2 = Targeter(xv2, fx2, maxiter, tol)
Xhist2, err2 = target(targ2)

println("After retargeting")
println(norm(tofullvector(fx2)))

println("\n9:2 NRHO")
nrho92_traj = Trajectory(model, [X1_fix, X2, X3], [T1, T2, T3])
println("Periodic? $(isperiodic(nrho92_traj))")
println("Continuous? $(iscontinuous(nrho92_traj))")

update!(X1, tofullvector(X1_fix))
