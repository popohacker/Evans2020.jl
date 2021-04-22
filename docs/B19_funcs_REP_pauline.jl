#YOOOO This one is Circe's branch

using Distributions
using LinearAlgebra

#questions 
# 1 : args python vs julia
# 2 : np. python ?? is inf... 

function trunc_norm_draws(unif_vals::Array, mu::Float64, sigma::Float64, cut_lb::Any, cut_ub::Any)
    if cut_lb == nothing && cut_ub == nothing
        cut_ub_cdf = 1.0
        cut_lb_cdf = 0.0
    elseif cut_lb != nothing && cut_ub == nothing
        cut_ub_cdf = 1.0
        cut_lb_cdf = cdf.(Normal(mu, sigma),  cut_lb)
    elseif cut_lb == nothing && cut_ub != nothing
        cut_ub_cdf = cdf.(Normal(mu, sigma), cut_ub)
        cut_lb_cdf = 0.0
    elseif cut_lb != nothing && cut_ub != nothing
        cut_ub_cdf = cdf.(Normal(mu, sigma), cut_ub)
        cut_lb_cdf = cdf.(Normal(mu, sigma),  cut_lb)
    end
	
    unif2_vals = transpose(reduce(hcat, unif_vals)) .* (cut_ub_cdf - cut_lb_cdf) .+ cut_lb_cdf # CHECK LATER -- on utilisait une matrice au format python pour tester la fonction mais peut être que quand elle sera directement générée par Julia pas besoin de transpose(reduce())
    tnorm_draws = quantile.(Normal(mu, sigma), unif2_vals)

    return tnorm_draws
end

function get_Y(k2t, zt, args)
    nvec, epsilon, alpha = args
    close_tol = 1e-6
    Kt = k2t
    Lt = sum(nvec) 
    At = exp.(zt)
    if isapprox(epsilon, 1.0, atol=close_tol) == true
        Yt = At * ((Kt)^alpha)*((Lt)^(1-alpha))
    elseif isinf(epsilon) == true
        Yt = At * (alpha * Kt + (1 - alpha) * Lt)
    elseif epsilon > 0 & isapprox(epsilon, 1.0, atol=close_tol) == false & isinf(epsilon) == false
        Yt = At * (alpha*(Kt^((epsilon -1)/epsilon)) + (1-alpha) *(Lt^((epsilon-1)/epsilon)))^(epsilon/(epsilon-1))
    elseif epsilon <= 0 
        err_msg = "ERROR get_Y(): epsilon <=0"
        print(err_msg)
    end

return Yt
end

function get_C(c1t, c2t)
    C = c1t + c2t
    return C  
end

function get_w(k2t, zt, args)
    nvec, epsilon, alpha = args
    Lt = sum(nvec)
    At = exp.(zt)
    if isinf(epsilon) == true
        wt = (1-alpha)*At
    elseif epsilon > 0 & isinf(epsilon) == false
        Yt = get_Y(k2t, zt, args)
        wt = (1-alpha)*(At^((epsilon-1)/epsilon))*((Yt/Lt)^(1/epsilon))
    elseif epsilon <=0
        err_msg = "ERROR get_w(): epsilon <=0"
        print(err_msg)
    end
return wt
end    

function get_r(k2t, zt, args)
    nvec, epsilon, alpha, delta = args
    Kt = k2t
    At = exp.(zt)
    if isinf(epsilon) == true
        rt = alpha * At - delta
    elseif epsilon > 0 & isinf(epsilon) == false
        Y_args = (nvec, epsilon, alpha)
        Yt = get_Y(k2t, zt, Y_args)
        rt = alpha * At^((epsilon-1)/epsilon)*(Yt/Kt)^(1/epsilon) - delta
    elseif epsilon <= 0 
        err_msg = "ERROR get_r(): epsilon <= 0"
        print(err_msg)
    end
return rt
end



