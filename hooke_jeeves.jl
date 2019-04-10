#! algorithm 7.6 of Kochenderfer & Wheeler (2019)

# α starting step size

basis(i, n) = [k == i ? 1.0 : 0.0 for k in 1 : n]

function hooke_jeeves(f, x, α, ϵ, γ=0.5)
    y, n = f(x), length(x)
    while α > ϵ
        improved = false
        x_best, y_best = x, y
        for i in 1:n
        	for sgn in (-1,1)
	            x′ = x + sgn*α*basis(i, n)
	            y′ = f(x′)
	            if y′ < y_best
	                x_best, y_best, improved = x′, y′, true
	            end
            end
        end
        x, y = x_best, y_best
        if !improved
            α *= γ
        end
    end
    return x
end
