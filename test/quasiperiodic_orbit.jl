using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test
using BenchmarkTools

######## Generate a QPO by targeting Jacobi constant ########
include("nrho92.jl");
ref = nrho92

# Longitudinal angle of PO at which to target the QPO
thT = pi

# Number of nodes on IC
N = 35

# Stepoff distance
d = 500

# Xstar
fixedpt = ref(thT)

############### FREE VARIABLE SETUP ############### 
# Generate FreeVariable for invariant curve
u0vec = linear_invariant_curve_2d(ref, thT, N, d)[1]

rmind = [1]
U0 = FreeVariable("U0", u0vec, includeinds=setdiff(1:length(u0vec),rmind))

# Generate FreeVariable for stroboscopic time
T = FreeVariable("T", period(ref))

# Generate FreeVariable for twist angle
lam, vee = center_eigs(ref, thT)
rho = FreeVariable("rho", angle(lam[1]))

# Generate XVector
xv = XVector(U0, T, rho)

############### CONSTRAINT SETUP ############### 
# Invariance Constraint
ic = InvarianceConstraint2D(U0, fixedpt, T, rho, model)

# Jacobi Constraint
jcc = JacobiConstraint(model, jacobi_constant(ref), U0; refpt=fixedpt);

# Generate FXVector
fx = FXVector(ic, jcc)

############### TARGETER SETUP ############### 
maxiter = 25
tol = 1e-10
targ = Targeter(xv, fx, maxiter, tol);

# Target QPO
Xhist, err = target(targ,debug=true);

############### QUASIPERIODIC ORBIT OBJECT ############### 
u0vec = tofullvector(U0)
q0 = similar(u0vec)
for i = 1:N
    q0[6*i-5:6*i] = u0vec[6i-5:6i] + fixedpt
end

qpo92 = QuasiPeriodicOrbit(model, q0, tofullvector(T)[1], 
                           tofullvector(rho)[1], fixedpt, tol=tol)
