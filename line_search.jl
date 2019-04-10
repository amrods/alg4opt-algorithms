#! algorithm 4.1 of Kochenderfer & Wheeler (2019)

function line_search(f, x, d)
    if norm(d) ≈ 0; return x; end; objective = α -> f(x + α*d)
    a, b = bracket_minimum(objective)
    α = minimize(objective, a, b)
    return x + α*d
end
