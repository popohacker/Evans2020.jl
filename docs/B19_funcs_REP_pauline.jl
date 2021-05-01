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

#### FROM HERE TO GET1PRMUC2PDF -- NO SYNTAX PROBLEM
#### I HAVENT CHECKED THAT IT GIVES THE SAMES RESULTS AS PYTHON YET

function get_Hbar_err(zt, x...)
    #
    # This function is the error function that solves for the current
    # period shock that sets w * n1 + x1 - c_min - K_min = Hbar. This is
    # the minimum shock that does not create default.
    # '''
    
    k2t, nvec, epsilon, alpha, Hbar, x1, c_min, K_min = x
    n1 = nvec[1]
    w_args = (nvec, epsilon, alpha)
    wt = get_w(k2t, zt, w_args)
    Hbar_err = Hbar - wt * n1 - x1 + c_min + K_min

    return Hbar_err
end



using Roots
using SymEngine

function get_zstar(k2t, ztm1, args)
    mu, rho, nvec, epsilon, alpha, Hbar, x1, c_min, K_min, sigma = args
    #mu, rho, nvec, epsilon, alpha, Hbar, x1, c_min, K_min, sigma = args
    z_init = 1.5 * mu
    z_mu = rho * ztm1 + (1 - rho) * mu
    zst_args = (k2t, nvec, epsilon, alpha, Hbar, x1, c_min, K_min)
    @vars zst_args
    results = find_zero(get_Hbar_err, z_init)
    z_star = results[1]
    eps_star = z_star - rho * ztm1 - (1 - rho) * mu
    A_star = exp.(z_star)
    prob_shut = cdf.(Normal(z_mu, sigma), z_star)
    if results.convergence_failed  
        err_msg = ("zstar ERROR: Root finder did not solve in ... get_zstar().")
        print("z_star = $z_star")
        print("Hbar_err = $(get_Hbar_err[1])")
        print(err_msg)
    end
    return z_star, eps_star, A_star, prob_shut
end


function get_c2t(k2t, zt, args)
    nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min = args
    w_args = (nvec, epsilon, alpha) 
    wt = get_w(k2t, zt, w_args)
    r_args = (nvec, epsilon, alpha, delta)
    rt = get_r(k2t, zt, r_args)
    H_args = (tau, Hbar, n1, x1, c_min, K_min)
    Ht, default = get_Ht(wt, H_args)
    c2t = (1 + rt) * k2t + wt * n2 + Ht

    return c2t
end

function get_MUc_CRRA(c, gamma)
    ###
    #--------------------------------------------------------------------
    #Generate marginal utility(ies) of consumption with CRRA consumption
    #utility and stitched function at lower bound such that the new
    #hybrid function is defined over all consumption on the real
    #line but the function has similar properties to the Inada condition.

    #u'(c) = c ** (-sigma) if c >= epsilon
    #      = g'(c) = 2 * b2 * c + b1 if c < epsilon

    #    such that g'(epsilon) = u'(epsilon)
    #    and g''(epsilon) = u''(epsilon)

    #    u(c) = (c ** (1 - sigma) - 1) / (1 - sigma)
    #    g(c) = b2 * (c ** 2) + b1 * c + b0
    #--------------------------------------------------------------------
    #INPUTS:
    #c  = scalar, individual consumption in a particular period
    #gamma = scalar >= 1, coefficient of relative risk aversion for CRRA
    #        utility function: (c**(1-gamma) - 1) / (1 - gamma)
    #graph = boolean, =True if want plot of stitched marginal utility of
    #        consumption function

    #OTHER FUNCTIONS AND FILES CALLED BY THIS FUNCTION: None

    #OBJECTS CREATED WITHIN FUNCTION:
    #epsilon    = scalar > 0, positive value close to zero
    #c_s        = scalar, individual consumption
    #c_s_cnstr  = boolean, =True if c_s < epsilon
    #b1         = scalar, intercept value in linear marginal utility
    #b2         = scalar, slope coefficient in linear marginal utility
    #MU_c       = scalar or (p,) vector, marginal utility of consumption
    #             or vector of marginal utilities of consumption
    #p          = integer >= 1, number of periods remaining in lifetime
    #cvec_cnstr = (p,) boolean vector, =True for values of cvec < epsilon

    #FILES CREATED BY THIS FUNCTION:
    #    MU_c_stitched.png

    #RETURNS: MU_c
    #--------------------------------------------------------------------
    ###
    c_epsilon = 1e-5
    if c > c_epsilon
        MUc = c .* (-gamma)
    elseif c <= c_epsilon
        b2 = (-gamma * (c_epsilon .* (-gamma - 1))) / 2
        b1 = (c_epsilon .* (-gamma)) - 2 * b2 * c_epsilon
        MUc = 2 * b2 * c + b1
    end
        # print('c=', c, ', MUc=', MUc)
    return MUc
end

function get_c1mgam(c, gamma)
    ###
    #--------------------------------------------------------------------
    #Generate marginal utility(ies) of consumption with CRRA consumption
    #utility and stitched function at lower bound such that the new
    #hybrid function is defined over all consumption on the real
    #line but the function has similar properties to the Inada condition.

    #f(c) = c ** (1-sigma) if c >= epsilon
    #g(c) = b2 * c + b1    if c < epsilon

    #    such that g(epsilon) = f(epsilon)
    #    and g'(epsilon) = f'(epsilon)

    #    f(c) = c ** (1 - sigma)
    #    g(c) = b2 * c + b1

    #    s.t. b2 = (1 - gamma) * (epsilon ** (-gamma))
    #         b1 = epsilon**(-gamma) - (1-gamma) * (epsilon ** (1-gamma))
    #--------------------------------------------------------------------
    #INPUTS:
    #c  = scalar, individual consumption in a particular period
    #gamma = scalar >= 1, coefficient of relative risk aversion for CRRA
    #        utility function: (c**(1-gamma) - 1) / (1 - gamma)
    #graph = boolean, =True if want plot of stitched marginal utility of
    #        consumption function

    #OTHER FUNCTIONS AND FILES CALLED BY THIS FUNCTION: None

    #OBJECTS CREATED WITHIN FUNCTION:
    #epsilon    = scalar > 0, positive value close to zero
    #b1         = scalar, intercept value in linear marginal utility
    #b2         = scalar, slope coefficient in linear marginal utility
    #MU_c       = scalar or (p,) vector, marginal utility of consumption
    #             or vector of marginal utilities of consumption
    #p          = integer >= 1, number of periods remaining in lifetime
    #cvec_cnstr = (p,) boolean vector, =True for values of cvec < epsilon

    #FILES CREATED BY THIS FUNCTION:
    #    MU_c_stitched.png

    #RETURNS: f_c
    #--------------------------------------------------------------------
    ###
    c_epsilon = 1e-5
    if c > c_epsilon
        f_c = c .* (1 - gamma)
    elseif c <= c_epsilon
        b2 = (1 - gamma) * (c_epsilon .* (-gamma))
        b1 = (c_epsilon .* (-gamma)) - b2 * c_epsilon
        f_c = b2 * c + b1
    end
    # print('c=', c, ', MUc=', MUc)
    return f_c
end

# @numba.jit(forceobj=True)
#he does that for paralellisation purposes

function LN_pdf(xvals, mu, sigma)
    ###
    #--------------------------------------------------------------------
    #This function gives the PDF of the lognormal distribution for xvals
    #given mu and sigma

    #(LN): f(x; mu, sigma) = (1 / (x * sigma * sqrt(2 * pi))) *
    #        exp((-1 / 2) * (((log(x) - mu) / sigma) ** 2))
    #        x in [0, infty), mu in (-infty, infty), sigma > 0
    #--------------------------------------------------------------------
    #INPUTS:
    #xvals = (N,) vector, data
    #mu    = scalar, mean of the ln(x)
    #sigma = scalar > 0, standard deviation of ln(x)

    #OTHER FUNCTIONS AND FILES CALLED BY THIS FUNCTION: None

    #OBJECTS CREATED WITHIN FUNCTION:
    #pdf_vals        = (N,) vector, probability of each observation given
    #                  the parameter values

    #FILES CREATED BY THIS FUNCTION: None

    #RETURNS: pdf_vals
    #--------------------------------------------------------------------
    ###
    pdf_vals = Array{Float64}(((1 / (sqrt(2 * π) * sigma * xvals)) *
                    exp((-1.0 / 2.0) * (((log.(xvals) - mu) / sigma) .* 2))))

    return pdf_vals
end

#@numba.jit(forceobj=True) 
# same here

function get_1pr_MU_c2_pdf(Atp1, x...)
    ###
    #This function is the target for calculating the integral
    #(expectation): E[(1+r_{tp1})*(c_{2,t+1})**(-gamma)]. This function
    #returns the value of
    #(1 + r_{tp1})*((c_{2,t+1})**(-gamma)) * pdf(A|mu,sigma)
    #for a given value of A and k2tp1
    ###
    (k2tp1, zt, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta, tau,
     Hbar, x1, c_min, K_min, gamma, sigma) = x
    ztp1 = log.(Atp1)
    z_mu = rho * zt + (1 - rho) * mu
    c2_args = (nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min)
    c2tp1 = get_c2t(k2tp1, ztp1, c2_args)
    r_args = (nvec, epsilon, alpha, delta)
    rtp1 = get_r(k2tp1, ztp1, r_args)
    MU_CRRA_c2tp1 = get_MUc_CRRA(c2tp1, gamma)
    MU_c2tp1_pdf = ((1 + rtp1) * MU_CRRA_c2tp1 *
                    (LN_pdf(Atp1, z_mu, sigma) / (1 - A_min_cdf)))
    return MU_c2tp1_pdf
end


#@numba.jit(forceobj=True)

function get_MU_c2_pdf(Atp1, x...)
    ###
    #This function is the target for calculating the integral
    #(expectation): E[(c_{2,t+1})**(-gamma)]. This function returns the
    #value of ((c_{2,t+1})**(-gamma)) * pdf(A|mu,sigma)
    #for a given value of A and k2tp1
    #'''
    (k2tp1, zt, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta, tau,
    Hbar, x1, c_min, K_min, gamma, sigma) = x
    ztp1 = log.(Atp1)
    z_mu = rho * zt + (1 - rho) * mu
    c2_args = (nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min)
    c2tp1 = get_c2t(k2tp1, ztp1, c2_args)
    MU_CRRA_c2tp1 = get_MUc_CRRA(c2tp1, gamma)
    MU_c2tp1_pdf = (MU_CRRA_c2tp1 *
                    (LN_pdf(Atp1, z_mu, sigma) / (1 - A_min_cdf)))

    return MU_c2tp1_pdf
end

#@numba.jit(forceobj=True)

function get_c2tp1_1mgam_pdf(Atp1, x...)
    ###
    #This function is the target for calculating the integral
    #(expectation): E[(c_{2,t+1})**(1-gamma)]. This function returns the
    #value of ((c_{2,t+1})**(1-gamma)) * pdf(A|mu,sigma)
    #for a given value of A and k2tp1
    ###
    (k2tp1, zt, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta, tau,
    Hbar, x1, c_min, K_min, gamma, sigma) = x
   ztp1 = log.(Atp1)
   z_mu = rho * zt + (1 - rho) * mu
   c2_args = (nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min)
   c2tp1 = get_c2t(k2tp1, ztp1, c2_args)
   c2tp1_1mgam = get_c1mgam(c2tp1, gamma)
   c2tp1_1mgam_pdf = (c2tp1_1mgam *
                      (LN_pdf(Atp1, z_mu, sigma) / (1 - A_min_cdf)))

   return c2tp1_1mgam_pdf
end

using QuadGK

function get_ExpMU_c2_b(k2pt1, zt, args)
    (k2tp1, zt, A_min, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta, tau,
     Hbar, x1, c_min, K_min, gamma, sigma) = args
     Ex_args = (k2tp1, zt, A_min_cdf, rho, mu, nvec, epsilon, alpha,
     delta, tau, Hbar, x1, c_min, K_min, gamma, sigma)
    (Exp_MU_CRRA_c2, _) = quadgk(Ex_args -> get_MU_c2_pdf(Ex_args), A_min, Inf)
    (Exp_c2_1mgam, _) = quandgk(Ex_args -> get_c2tp1_1mgam_pdf(Ex_args), A_min, Inf)
    MU = Exp_MU_CRRA_c2 / Exp_c2_1mgam
return MU
end

function get_Eul_err(k2tp1, x...)
    (k2t, zt, Ht, nvec, epsilon, beta, alpha, delta, x1, rho, mu, sigma,
     A_min, tau, Hbar, c_min, K_min, gamma, sigma) = x 
    n1 = nvec[1]
    w_args = (nvec, epsilon, alpha)
    wt = get_w(k2t, zt, w_args)
    c1 = wt * n1 + x1 - k2tp1 - Ht
    MU_c1 = get_MUc_CRRA(c1, 1.0)
    mu_ztp1 = rho * zt + (1 - rho) * mu
    if A_min == 0.0
        A_min_cdf = 0.0
    elseif A_min > 0.0
        A_min_cdf = cdf.(Normal(mu_ztp1, sigma), log(A_min))
    end
    MU_args = (A_min, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta,
    tau, Hbar, x1, c_min, K_min, gamma, sigma)
    Exp_MU_ctp1 = get_ExpMU_c2tp1_k(k2tp1, zt, MU_args)
    Eul_err = MU_c1 - (beta / (1 - beta)) * Exp_MU_ctp1
    return Eul_err
end


