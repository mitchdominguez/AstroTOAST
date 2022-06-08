using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test

# Generate a QPO by targeting Jacobi constant

# Obtain the 2:1 Halo
include("halo21.jl");
ref = halo21

# include("nrho92.jl");
# ref = nrho92

# Longitudinal angle of PO at which to target the QPO
thT = pi

# Number of nodes on IC
N = 35

# Stepoff distance
d = 300

# Xstar
fixedpt = ref(thT)

############### FREE VARIABLE SETUP ############### 
# Generate FreeVariable for invariant curve
u0vec = linear_invariant_curve_2d(ref, thT, N, d)[1]

rmind = []
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

Xhist, err = target(targ);

err
