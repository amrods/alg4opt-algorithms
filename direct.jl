#! algorithm 7.8 of Kochenderfer & Wheeler (2019)

function direct(f, a, b, ϵ, k_max)
    g = reparameterize_to_unit_hypercube(f, a, b)
    intervals = Intervals()
    n = length(a)
    c = fill(0.5, n)
    interval = Interval(c, g(c), fill(0, n))
    add_interval!(intervals, interval)
    c_best, y_best = copy(interval.c), interval.y

    for k in 1:k_max
        S = get_opt_intervals(intervals, ϵ, y_best)
        to_add = Interval[]
        for interval in S
            append!(to_add, divide(g, interval))
            dequeue!(intervals[min_depth(interval)])
        end
        for interval in to_add
            add_interval!(intervals, interval)
            if interval.y < y_best
                c_best, y_best = copy(interval.c), interval.y
            end
        end
    end

    return rev_unit_hypercube_parameterization(c_best, a, b)
end

rev_unit_hypercube_parameterization(x, a, b) = x .* (b - a) + a
function reparameterize_to_unit_hypercube(f, a, b)
    Δ = b - a
    return x -> f(x.*Δ + a)
end

using DataStructures

struct Interval
    c
    y
    depths
end

min_depth(interval) = minimum(interval.depths)

const Intervals = Dict{Int, PriorityQueue{Interval, Float64}}

function add_interval!(intervals, interval)
	d = min_depth(interval)
    if !haskey(intervals, d)
        intervals[d] = PriorityQueue{Interval, Float64}()
    end

    return enqueue!(intervals[d], interval, interval.y)
end

function get_opt_intervals(intervals, ϵ, y_best)
    max_depth = maximum(keys(intervals))
    stack = [DataStructures.peek(intervals[max_depth])[1]]
    d = max_depth-1
    while d ≥ 0
        if haskey(intervals, d) && !isempty(intervals[d])
            interval = DataStructures.peek(intervals[d])[1]
            x, y = 0.5*3.0^(-min_depth(interval)), interval.y

            while !isempty(stack)
            	interval1 = stack[end]
            	x1 = 0.5*3.0^(-min_depth(interval1))
            	y1 = interval1.y
            	l1 = (y - y1)/(x - x1)
            	if y1 - l1*x1 > y_best - ϵ || y < y1
            		pop!(stack)
            	elseif length(stack) > 1
            		interval2 = stack[end-1]
            		x2 = 0.5*3.0^(-min_depth(interval2))
            		y2 = interval2.y
            		l2 = (y1 - y2)/(x1 - x2)
            		if l2 > l1
            			pop!(stack)
                    else
                        break
            		end
                else
                    break
            	end
            end

            push!(stack, interval) # add new point
        end
        d -= 1
    end
    return stack
end

basis(i, n) = [k == i ? 1.0 : 0.0 for k in 1 : n]

function divide(f, interval)
    c, d, n = interval.c, min_depth(interval), length(interval.c)
    dirs = findall(interval.depths .== d)
    cs = [(c + 3.0^(-d-1)*basis(i,n),
           c - 3.0^(-d-1)*basis(i,n)) for i in dirs]
    vs = [(f(C[1]), f(C[2])) for C in cs]
    minvals = [min(V[1], V[2]) for V in vs]

    intervals = Interval[]
    depths = copy(interval.depths)
    for j in sortperm(minvals)
        depths[dirs[j]] += 1
        C, V = cs[j], vs[j]
        push!(intervals, Interval(C[1], V[1], copy(depths)))
        push!(intervals, Interval(C[2], V[2], copy(depths)))
    end
    push!(intervals, Interval(c, interval.y, copy(depths)))
    return intervals
end
