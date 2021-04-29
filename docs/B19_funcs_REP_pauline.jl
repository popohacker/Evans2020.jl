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

function get_Ht(wt,args)
	tau, Hbar, n1, x1, c_min, K_min = args
	default = "false"
	if tau == nothing 
		Ht = min(Hbar, wt * n1 + x1 - c_min - K_min)
		if Ht < Hbar
			default = "True"
	elseif tau > 0 && tau < 1
			Ht = tau * (wt * n1 + x1)
	elseif tau <= 0 || tau >= 1
		err_msg = "ERROR get_Ht(): tau is not valid value."
        print(err_msg)
		end
	end
	return Ht, default	
end

#new changes - Pauline
function get_Hbar_err(zt, length(args))
    #
    # This function is the error function that solves for the current
    # period shock that sets w * n1 + x1 - c_min - K_min = Hbar. This is
    # the minimum shock that does not create default.
    # '''
    k2t, nvec, epsilon, alpha, Hbar, x1, c_min, K_min = args
    n1 = nvec[1]
    w_args = (nvec, epsilon, alpha)
    wt = get_w(k2t, zt, w_args)
    Hbar_err = Hbar - wt * n1 - x1 + c_min + K_min

    return Hbar_err
end

using Roots

function get_zstar(k2t, ztm1, args)
    mu, rho, nvec, epsilon, alpha, Hbar, x1, c_min, K_min, sigma = args
    mu, rho, nvec, epsilon, alpha, Hbar, x1, c_min, K_min, sigma = args
    z_init = 1.5 * mu
    z_mu = rho * ztm1 + (1 - rho) * mu
    zst_args = (k2t, nvec, epsilon, alpha, Hbar, x1, c_min, K_min)
    @vars zst_args
    results = find_zero(get_Hbar_err, z_init)
    z_star = results[1]
    eps_star = z_star - rho * ztm1 - (1 - rho) * mu
    A_star = exp.(z_star)
    prob_shut = cdf.(Normal(z_mu, sigma), z_star)
    if results.convergence_failed = true 
        err_msg = ("zstar ERROR: Root finder did not solve in" +
        "get_zstar().")
        print("z_star = $z_star")
        print("Hbar_err = $(get_Hbar_err[1])")
        print(err_msg)
    end
    return z_star, eps_star, A_star, prob_shut
end
