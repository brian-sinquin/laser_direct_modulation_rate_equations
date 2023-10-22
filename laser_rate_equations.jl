using DifferentialEquations

function laser!(du, u, p, s)

    (r, m, η, ε, α) = p # parameters

    e = u[1] # optical field
    n = u[2] # carrier inversion

    Jm = r * (1 + m * cos(2π * η * s))
    J = Jm * (Jm >= 0) # clamping  # DC + AC pump must be null or positive (a laser is a diode)

    du[1] = 0.5(1 + 1im * α) * n * e # E. field evolution
    du[2] = ε * (J - (n + 1) * (1 + abs2(e))) # Inversion evolution

end


function integrate(p, ta, tb, Np)

    r, m, f0, τp, τc, α = p

    η = f0 * τp # normalized frequency
    ε = τp / τc # lifetime ratio

    p0 = (r, m, η, ε, α)

    X0 = [0.000001 + 0.000001im, 0.000001] # initial condition 
    s = LinRange(ta, tb, Np) ./ τp
    prob = ODEProblem(laser!, X0, (0.0, s[end]), p0)
    sol = solve(prob, RK4(), saveat=(tb - ta) / τp / Np)
    sol = sol(s)
    t = s .* τp
    E = sol[1, :]
    n = real.(sol[2, :])

    return E, n, t

end