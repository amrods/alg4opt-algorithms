#! algorithm 4.4 of Kochenderfer & Wheeler (2019)

function trust_region_descent(f, ∇f, H, x, k_max; η1=0.25, η2=0.5, γ1=0.5, γ2=2.0, δ=1.0)
	y = f(x)
	for k in 1 : k_max
		x′, y′ = solve_trust_region_subproblem(∇f, H, x, δ)
		r = (y - f(x′)) / (y - y′)
		if r < η1
			δ *= γ1
		else
			x, y = x′, y′
			if r > η2
				δ *= γ2
			end
		end
	end
	return x
end

using Convex
using SCS

const ⋅ = dot

function solve_trust_region_subproblem(∇f, H, x0, δ)
	x = Variable(length(x0))
	p = Convex.minimize(∇f(x0)⋅(x-x0) + quadform(x-x0, H(x0))/2)
	p.constraints += norm(x-x0) <= δ
	solve!(p, SCSSolver(verbose=false), verbose=false)
	return (x.value, p.optval)
end
