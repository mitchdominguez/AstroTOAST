using AstroTOAST
using LinearAlgebra
using StaticArrays
using Test
using BenchmarkTools

######## Generate a QPO by targeting Jacobi constant ########
# Obtain the 2:1 Halo
# include("halo21.jl");
# ref = halo21

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
@time Xhist, err = target(targ,debug=true);

# dim_state_0 = dimensionalize_state(dimensional_quantity_set(model), U0[1:6])
# println("\nDimensional State 0: $(dim_state_0)")
#
u0vec = tofullvector(U0)
q0 = similar(u0vec)
for i = 1:N
    q0[6*i-5:6*i] = u0vec[6i-5:6i] + fixedpt
end

qpo92 = QuasiPeriodicOrbit(model, q0, tofullvector(T)[1], 
                           tofullvector(rho)[1], fixedpt, tol=tol;
                           name="quasi92", family="9:2 Quasihalos", 
                           thT_offset=time2angle(tofullvector(T)[1], angle2time(ref, thT)))




@testset "qposs.jl" begin
    @test length(err) == 10
    @test norm(tofullvector(fx)) < tol
    @test norm(fx) ≈ 2.869185571417227e-11 atol=1e-12
    @test typeof(qpo92) == QuasiPeriodicOrbit{6,35}
end
