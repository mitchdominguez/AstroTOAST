#=
# lagrange.jl
#
# Capabilities for calculating lagrange point locations.
#
# Rolfe Power
=#

module LagrangePointHelpers

"""
    lp_collinear_approx_x(mu)

Calculate approximate `x` locations for collinear Lagrange points
"""
function lp_collinear_approx_x(mu)
    eta = (mu / 3.0)^(1.0/3.0);
    sigma = 7.0 * mu/12.0;

    xl1 = 1.0 - mu -
        eta * (1.0 - eta/3.0 - eta^2/9.0 -
            23.0 * eta^3/81.0 + 151.0 * eta^4/243.0 - eta^5/9.0)
    xl2 = 1.0 - mu +
        eta * (1.0 + eta/3.0 - eta^2/9.0 -
            31.0 * eta^3/81.0 - 119.0 * eta^4/243.0 - eta^5/9.0)
    xl3 = -mu - 1.0 +
        sigma * (1.0 + 23.0 * sigma^2/84.0 +
            23.0 * sigma^3/84.0 + 761.0 * sigma^4/2352.0 +
            3163.0 * sigma^5/7056.0 + 30703.0 * sigma^6/49392.0)

    (xl1, xl2, xl3)
end

"""
    eq1(mu, g)

Helper function for running Newton-Raphson on L1 Lagrange point.
"""
function eq1(mu, g)
    f = 1-mu-g-(1-mu)/(1-g)^2 + mu/g^2
    df = -1 - 2*(1-mu)/(1-g)^3 - 2*mu/g^3
    (f, df)
end

"""
    eq2(mu, g)

Helper function for running Newton-Raphson on L2 Lagrange point.
"""
function eq2(mu, g)
    f = (1-mu)/(1+g)^2 + mu/g^2 - 1 + mu - g
    df = -2*(1-mu)/((1+g)^3) - 2*mu/g^3 - 1
    (f, df)
end

"""
    eq3(mu, g)

Helper function for running Newton-Raphson on L3 Lagrange point.
"""
function eq3(mu, g)
    f = -mu - g + (1-mu)/g^2 + mu/(g+1)^2
    df = -1 - 2*(1-mu)/g^3 - 2*mu/(g+1)^3
    (f, df)
end

end # module LagrangePointHelpers

"""
    equilibrium_solutions(m::Cr3bpModel; [max_iter=100], [tolerance=1e-12])

Calculate the Lagrange point locations for the CR3BP
"""
function equilibrium_solutions(m::Cr3bpModel; max_iter=100, tolerance=1e-12)
    mu = mass_ratio(m)

    f1(g) = LagrangePointHelpers.eq1(mu, g)
    f2(g) = LagrangePointHelpers.eq2(mu, g)
    f3(g) = LagrangePointHelpers.eq3(mu, g)

    (x10, x20, x30) = LagrangePointHelpers.lp_collinear_approx_x(mu)

    g1 = 1.0 - mu - x10
    g2 = x20 - 1.0 + mu
    g3 = -mu - x30

    # x1
    count = 0
    f, df = f1(g1)
    while abs(f) > tolerance && count < max_iter
        g1 -= f / df
        f, df = f1(g1)
        count += 1
    end
    x1 = if count < max_iter
        1.0 - mu - g1
    else
        NaN
    end

    # x2
    count = 0
    f, df = f2(g2)
    while abs(f) > tolerance && count < max_iter
        g2 -= f / df
        f, df = f2(g2)
        count += 2
    end
    x2 = if count < max_iter
        1.0 - mu + g2
    else
        NaN
    end

    # x3
    count = 0
    f, df = f3(g3)
    while abs(f) > tolerance && count < max_iter
        g3 -= f / df
        f, df = f3(g3)
        count += 3
    end
    x3 = if count < max_iter
        -mu - g3
    else
        NaN
    end

    x45 = 0.5 - mu;
    y45 = sqrt(3)/2

    L1 = @SVector [  x1, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    L2 = @SVector [  x2, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    L3 = @SVector [  x3, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    L4 = @SVector [ x45, y45, 0.0, 0.0, 0.0, 0.0 ]
    L5 = @SVector [ x45,-y45, 0.0, 0.0, 0.0, 0.0 ]

    @SVector [L1, L2, L3, L4, L5]
end

"""
    lyapunov_linear_ics(model::Cr3bpModel, xi_0::Real, eta_0::Real, index::Int)

Calculate the linear approximation for the Lyapunov orbit around the 
collinear lagrange points. `index` determines which collinear lagrange point
the initial conditions will be taken about, and so it must be equal to either
1,2, or 3.
"""
function lyapunov_linear_ics(model::Cr3bpModel, xi_0::Real, eta_0::Real, index::Int)
    if index!=1 && index!=2 && index!=3
        throw(InvalidStateException("index must be equal to 1, 2, or 3", :index))
    end

    Li = equilibrium_solutions(model)[index]

    r_L = Li[1:3]

    r_0 = r_L + [xi_0,eta_0,0]

    U = pseudopotential_jacobian(model, Li)
    U_xx = U[1,1]
    U_yy = U[2,2]

    # Solve for beta_1, beta_2
    beta_1 = 2 - ((U_xx+U_yy)/2)
    beta_2 = sqrt(complex(-U_xx*U_yy))

    # Solve for s, beta_3
    s = sqrt(beta_1 + sqrt(beta_1^2 + beta_2^2))
    beta_3 = (s^2 + U_xx)/(2*s)

    # Position is r_L1 + [xi_0,eta_0]
    # Velocity is [xi_0_dot, eta_0_dot]

    xi_0_dot = eta_0*s/beta_3
    eta_0_dot = -beta_3*xi_0*s
    v_0 = [xi_0_dot,eta_0_dot,0]

    return IC = ignore_imag(vcat(r_0, v_0))
end
