using Distributions
using LinearAlgebra
using Roots
using SymEngine
using QuadGK
using Roots
using DifferentialEquations
using CommonSolve
using Optim
using LeastSquaresOptim 

#Function 1
function trunc_norm_draws(unif_vals::Any, mu::Any, sigma::Any, cut_lb::Any = nothing, cut_ub::Any = nothing)
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
    unif2_vals = transpose(reduce(hcat, unif_vals)) .* (cut_ub_cdf - cut_lb_cdf) .+ cut_lb_cdf 
    tnorm_draws = quantile.(Normal(mu, sigma), unif2_vals)
    return tnorm_draws
end

#Function 2
function get_Y(k2t, zt, args)
    nvec, epsilon, alpha = args
    close_tol = 1e-6
    Kt = k2t
    Lt = sum(nvec) 
    At = exp.(zt)
    if isapprox(epsilon, 1.0, atol=close_tol) == true
        Yt = At .* ((Kt).^alpha)*((Lt).^(1- alpha))
    elseif isinf(epsilon) == true
        Yt = At .* (alpha .* Kt .+ (1 - alpha) .* Lt)
    elseif epsilon > 0 & isapprox(epsilon, 1.0, atol=close_tol) == false & isinf(epsilon) == false
        Yt = At .* (alpha.*(Kt.^((epsilon -1)/epsilon)) + (1-alpha) *(Lt.^((epsilon-1)/epsilon))).^(epsilon/(epsilon-1))
    elseif epsilon <= 0 
        err_msg = "ERROR get_Y(): epsilon <=0"
        print(err_msg)
    end
return Yt
end

#Function 3
function get_C(c1t, c2t)
    C = c1t + c2t
    return C  
end

#Function 4
function get_w(k2t, zt, args)
    nvec, epsilon, alpha = args
    Lt = sum(nvec)
    At = exp.(zt)
    if isinf(epsilon) == true
        wt = (1-alpha).*At
    elseif epsilon > 0 & isinf(epsilon) == false
        Yt = get_Y(k2t, zt, args)
        wt = (1-alpha)*(At.^((epsilon-1)/epsilon)).*((Yt/Lt).^(1/epsilon))
    elseif epsilon <=0
        err_msg = "ERROR get_w(): epsilon <=0"
        print(err_msg)
    end
return wt
end    

#Function 5
function get_r(k2t, zt, args)
    nvec, epsilon, alpha, delta = args
    Kt = k2t
    At = exp.(zt)
    if isinf(epsilon) == true
        rt = alpha .* At .- delta
    elseif epsilon > 0 & isinf(epsilon) == false
        Y_args = (nvec, epsilon, alpha)
        Yt = get_Y(k2t, zt, Y_args)
        rt = alpha .* At.^((epsilon-1)/epsilon).*(Yt/Kt).^(1/epsilon) .- delta
    elseif epsilon <= 0 
        err_msg = "ERROR get_r(): epsilon <= 0"
        print(err_msg)
    end
return rt
end

#Function 6
function get_Ht(wt, args)
	tau, Hbar, n1, x1, c_min, K_min = args
    default = false
    if tau === nothing 
        x = (wt .* n1 .+ x1 .- c_min .- K_min)[1]
        Ht = min(Hbar, x)
		if Ht < Hbar
            default = true
        end
	elseif tau > 0 && tau < 1
            Ht = tau .* (wt .* n1 .+ x1)
	elseif tau <= 0 || tau >= 1
		err_msg = "ERROR get_Ht(): tau is not valid value."
        print(err_msg)
	end
	return Ht, default	
end

#Function 7
function get_Hbar_err(zt, args...)
    # This function is the error function that solves for the current
    # period shock that sets w * n1 + x1 - c_min - K_min = Hbar. This is
    # the minimum shock that does not create default.
    k2t, nvec, epsilon, alpha, Hbar, x1, c_min, K_min = args
    n1 = nvec[1]
    w_args = (nvec, epsilon, alpha)
    wt = get_w(k2t, zt, w_args)
    Hbar_err = Hbar - wt * n1 - x1 + c_min + K_min
    return Hbar_err
end

#Function 8
function get_zstar(k2t, ztm1, args)
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
    if results.convergence_failed 
        err_msg = ("zstar ERROR: Root finder did not solve in" +
        "get_zstar().")
        print("z_star = $z_star")
        print("Hbar_err = $(get_Hbar_err[1])")
        print(err_msg)
    end
    return z_star, eps_star, A_star, prob_shut
end

#Function 9
function get_Hbar_err(zt, x...)
    # This function is the error function that solves for the current
    # period shock that sets w * n1 + x1 - c_min - K_min = Hbar. This is
    # the minimum shock that does not create default.
    k2t, nvec, epsilon, alpha, Hbar, x1, c_min, K_min = x
    n1 = nvec[1]
    w_args = (nvec, epsilon, alpha)
    wt = get_w(k2t, zt, w_args)
    Hbar_err = Hbar - wt * n1 - x1 + c_min + K_min
    return Hbar_err
end

#Function 10
function get_zstar(k2t, ztm1, args)
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
    if results.convergence_failed 
        err_msg = ("zstar ERROR: Root finder did not solve in ... get_zstar().")
        print("z_star = $z_star")
        print("Hbar_err = $(get_Hbar_err[1])")
        print(err_msg)
    end
    return z_star, eps_star, A_star, prob_shut
end

#Function 11
function get_c2t(k2t, zt, args)
    nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min = args
    n1, n2 = nvec
    w_args = (nvec, epsilon, alpha) 
    wt = get_w(k2t, zt, w_args)
    r_args = (nvec, epsilon, alpha, delta)
    rt = get_r(k2t, zt, r_args)
    H_args = (tau, Hbar, n1, x1, c_min, K_min)
    Ht, default = get_Ht(wt, H_args)
    c2t = (1 .+ rt) .* k2t .+ wt .* n2 .+ Ht
    return c2t
end

#Function 12
function get_MUc_CRRA(c, gamma)
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

    c_epsilon = 1e-5
    if c[1] > c_epsilon 
        MUc = c.^(-gamma)
    elseif c[1] <= c_epsilon 
        b2 = (-gamma * (c_epsilon^(-gamma - 1))) / 2
        b1 = (c_epsilon^(-gamma)) .- 2 .* b2 .* c_epsilon
        MUc = 2 .* b2 .* c .+ b1
    end
    #print("c = ", c, ", MUc = ", MUc)
    global MUc
    return MUc
end

#Function 13
function get_c1mgam(c, gamma)
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
    c_epsilon = 1e-5
    if c[1] > c_epsilon
        f_c = c.^(1 - gamma)
    elseif c[1] <= c_epsilon
        b2 = (1 - gamma) * (c_epsilon^(-gamma))
        b1 = (c_epsilon^(-gamma)) - b2 * c_epsilon
        f_c = b2 * c .+ b1
    end
    #print("c = ", c, ", MUc = ", MUc)
    global f_c
    return f_c
end

#Function 14
function LN_pdf(xvals, mu, sigma)
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
    pdf_vals = (((1 / (sqrt(2 * ??) * sigma .* xvals)) .*
                   exp((-1.0 / 2.0) .* (((log.(xvals) .- mu) ./ sigma)^2))))
    return pdf_vals
end

#Function 15
function get_1pr_MU_c2_pdf(Atp1, args)
    #This function is the target for calculating the integral
    #(expectation): E[(1+r_{tp1})*(c_{2,t+1})**(-gamma)]. This function
    #returns the value of
    #(1 + r_{tp1})*((c_{2,t+1})**(-gamma)) * pdf(A|mu,sigma)
    #for a given value of A and k2tp1
    (k2tp1, zt, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min, gamma, sigma) = args
    ztp1 = log(Atp1)
    z_mu = rho * zt + (1 - rho) * mu
    c2_args = (nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min)
    c2tp1 = get_c2t(k2tp1, ztp1, c2_args)
    r_args = (nvec, epsilon, alpha, delta)
    rtp1 = get_r(k2tp1, ztp1, r_args)
    MU_CRRA_c2tp1 = get_MUc_CRRA(c2tp1, gamma)
    MU_c2tp1_pdf = ((1 .+ rtp1) .* MU_CRRA_c2tp1 .* (LN_pdf(Atp1, z_mu, sigma) / (1 .- A_min_cdf)))
    return MU_c2tp1_pdf
end

#Function 16
function get_MU_c2_pdf(Atp1, args)
    #This function is the target for calculating the integral
    #(expectation): E[(c_{2,t+1})**(-gamma)]. This function returns the
    #value of ((c_{2,t+1})**(-gamma)) * pdf(A|mu,sigma)
    #for a given value of A and k2tp1
    k2tp1, zt, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min, gamma, sigma = args
    ztp1 = log.(Atp1)
    z_mu = rho * zt + (1 - rho) * mu
    c2_args = (nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min)
    c2tp1 = get_c2t(k2tp1, ztp1, c2_args)
    MU_CRRA_c2tp1 = get_MUc_CRRA(c2tp1, gamma)
    MU_c2tp1_pdf = (MU_CRRA_c2tp1 * (LN_pdf(Atp1, z_mu, sigma) / (1 - A_min_cdf)))
    return MU_c2tp1_pdf
end

#Function 17
function get_c2tp1_1mgam_pdf(Atp1, args)
    #This function is the target for calculating the integral
    #(expectation): E[(c_{2,t+1})**(1-gamma)]. This function returns the
    #value of ((c_{2,t+1})^(1-gamma)) * pdf(A|mu,sigma)
    #for a given value of A and k2tp1
    k2tp1, zt, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min, gamma, sigma = args
    ztp1 = log.(Atp1)
    z_mu = rho * zt + (1 - rho) * mu
    c2_args = (nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min)
    c2tp1 = get_c2t(k2tp1, ztp1, c2_args)
    c2tp1_1mgam = get_c1mgam(c2tp1, gamma)
    c2tp1_1mgam_pdf = (c2tp1_1mgam * (LN_pdf(Atp1, z_mu, sigma) / (1 - A_min_cdf)))
    return c2tp1_1mgam_pdf
end

#Function 18
function get_ExpMU_c2tp1_k(k2tp1, zt, args)
    (A_min, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta, tau, Hbar, x1,  c_min, K_min, gamma, sigma) = args
    Ex_args = (k2tp1, zt, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min, gamma, sigma)
    (Exp_1pr_MU_CRRA_c2, _) = quadgk(x -> get_1pr_MU_c2_pdf(x, Ex_args), A_min, Inf)
    (Exp_c2_1mgam, _) = quadgk(x -> get_c2tp1_1mgam_pdf(x, Ex_args), A_min, Inf)
    MU = Exp_1pr_MU_CRRA_c2 / Exp_c2_1mgam
    return MU
end

#Function 19
function get_ExpMU_c2_b(k2pt1, zt, args)
    (k2tp1, zt, A_min, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta, tau,
    Hbar, x1, c_min, K_min, gamma, sigma) = args
    Ex_args = (k2tp1, zt, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min, gamma, sigma)
    (Exp_MU_CRRA_c2, _) = quadgk(x -> get_MU_c2_pdf(x, Ex_args), A_min, Inf)
    (Exp_c2_1mgam, _) = quadgk(x -> get_c2tp1_1mgam_pdf(x, Ex_args), A_min, Inf)
    MU = Exp_MU_CRRA_c2 / Exp_c2_1mgam
return MU
end

#Function 20
function get_Eul_err(k2tp1, args)
    (k2t, zt, Ht, nvec, epsilon, beta, alpha, delta, x1, rho, mu, sigma, A_min, tau, Hbar, c_min, K_min, gamma, sigma) = args
    n1 = nvec[1]
    w_args = (nvec, epsilon, alpha)
    wt = get_w(k2t, zt, w_args)
    c1 = wt .* n1 .+ x1 .- k2tp1 .- Ht
    MU_c1 = get_MUc_CRRA(c1, 1.0)
    mu_ztp1 = rho * zt + (1 - rho) * mu
    if A_min == 0.0
        A_min_cdf = 0.0
    elseif A_min > 0.0
        A_min_cdf = cdf.(Normal(mu_ztp1, sigma), log(A_min))
    end
    MU_args = (A_min, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min, gamma, sigma)
    Exp_MU_ctp1 = get_ExpMU_c2tp1_k(k2tp1, zt, MU_args)
    Eul_err = MU_c1 - (beta / (1 - beta)) * Exp_MU_ctp1
    global MU_c1
    return real(Eul_err[1,1])
end

#Function 21
function get_k2tp1(k2t, zt, args)
    #--------------------------------------------------------------------
    #Solve for k2tp1
    #c1t + k2tp1 = wt * n1 - tau * w1 * n1
    #--------------------------------------------------------------------

    (Hbar, beta, gamma, k20, rho, mu, sigma, x1, nvec, epsilon, alpha, delta, tau, c_min, K_min, A_min, yrs_in_per) = args
    n1 = nvec[1]
    krange_tol = 0.01
    w_args = (nvec, epsilon, alpha)
    wt = get_w(k2t, zt, w_args)
    H_args = (tau, Hbar, n1, x1, c_min, K_min)
    Ht, default = get_Ht(wt, H_args)
    r_args = (nvec, epsilon, alpha, delta)
    rt = get_r(k2t, zt, r_args)
    c2_args = (nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min)
    c2t = get_c2t(k2t, zt, c2_args)

    # Compute A_min_cdf
    mu_ztp1 = rho * zt + (1 - rho) * mu
    if A_min == 0.0
        A_min_cdf = 0.0
    elseif A_min > 0.0
        A_min_cdf = cdf.(Normal(mu_ztp1, sigma), log(A_min))
    end

    if default == true
        print("Default: Ht < Hbar ==> wt * n1 + x1 - c_min - K_min < Hbar")
        k2tp1 = K_min
        c1t = c_min
        MU_c1 = get_MUc_CRRA(c1t, 1.0)
        Eul_err = 0.0
    elseif default == false
        # wt * n1 + x1 - c_min - K_min >= Hbar
        # Compute k2tp1_max
        k2tp1_max = wt * n1 + x1 - c_min - Ht
        if (k2tp1_max - K_min) < krange_tol && (k2tp1_max - K_min) >= 0.0
            print("Too small maximization range: k2tp1_max - K_min too small.")
            k2tp1 = 0.5 * K_min + 0.5 * k2tp1_max
            c1t = wt * n1 + x1 - k2tp1 - Ht
            MU_c1 = get_MUc_CRRA(c1t, 1.0)
            eul_args = (k2t, zt, Ht, nvec, epsilon, beta, alpha, delta,
                        x1, rho, mu, sigma, A_min, tau, Hbar, #remove A_min_cdf to correct an error
                        c_min, K_min, gamma, sigma)
            Eul_err = get_Eul_err(k2tp1, eul_args)
        elseif (k2tp1_max - K_min) < 0.0
            print("ERROR MSG : Problem in get_k2tp1() : k2tp1_max - K_min <= 0.")
        elseif (k2tp1_max - K_min) >= krange_tol
            k2_init = 0.5 * k2tp1_max + 0.5 * K_min
            eul_args = (k2t, zt, Ht, nvec, epsilon, beta, alpha, delta,
                        x1, rho, mu, sigma, A_min, tau, Hbar, c_min,
                        K_min, gamma, sigma)
            #print("K_min = ", K_min, ", k2tp1_max = ", k2tp1_max)
            results = optimize(k2tp1 -> get_Eul_err(k2tp1, eul_args), [k2_init], LevenbergMarquardt(), lower = fill(K_min, length(k2_init)))
            k2tp1 = results.minimizer
            c1t = wt .* n1 .+ x1 .- k2tp1 .- Ht
            Eul_err = [0.] 
            MU_c1 = get_MUc_CRRA(c1t, 1.0)

            if results.converged == false
                print("Results")
                print(results)
                print(" ")
                print("eul_args")
                print(eul_args)
                print(" ")
                print(" ")
                err_msg = "Root finder did not solve in get_Eul_err()."
                #err_msg = ("Minimization did not solve in get_neg_lf_util().")
                print(err_msg)
            end
        end
    end
   
    # Compute price of riskless one-period bond
    global MU_c1
    global k2tp1
    global c1t
    global Ht
    global c2t
    global wt
    global rt
    global rbar_t
    global rbar_t_an
    global default
    global Eul_err
    Ex_args = k2tp1, zt, A_min, A_min_cdf, rho, mu, nvec, epsilon, alpha, delta, tau, Hbar, x1, c_min, K_min, gamma, sigma
    Exp_MU_c2tp1 = get_ExpMU_c2_b(k2tp1, zt, Ex_args)
    pbar_t = (beta / (1 - beta)) .* (Exp_MU_c2tp1 ./ MU_c1) 
    rbar_t = (1 ./ pbar_t) .- 1 
    rbar_t_an = ((1 ./ pbar_t).^(1/yrs_in_per)) .- 1
    return k2tp1, c1t, Ht, c2t, wt, rt, rbar_t, rbar_t_an, default, Eul_err
end 

#Function 22
function sim_timepath(
    Hbar, beta, gamma, k20, sigma, x1, T, z0, z_min, rho, mu, nvec,
    epsilon, alpha, delta, tau, c_min, K_min, A_min, yrs_in_per,
    H_ind, risk_type_ind, risk_val_ind,
    avgrtp1_ind, avgrbart_ind, S_ind, zt_vec,
    rand_seed)

    #(k20, z0, nvec, x1, x2, c_min, K_min, Hbar, tau, beta, gamma, alpha, eps, delta, mu, rho, sigma, A_min, z_min, yrs_in_per, T) = args

    if H_ind === nothing
        H_ind = 0
    end 
    if risk_type_ind  === nothing
        risk_type_ind = 0
    end 
    if risk_val_ind  === nothing
        risk_val_ind = 0
    end 
    if avgrtp1_ind === nothing
        avgrtp1_ind = 0
    end 
    if avgrbart_ind  === nothing
        avgrbart_ind = 0
    end 
    if S_ind  === nothing
        S_ind = 0
    end 

    if zt_vec  === nothing
        if rand_seed  === nothing
            rand_seed = rand(1:1000)
        end 

        zt_vec = zeros(T)
        Random.seed!(rand_seed)
        unif_vec = rand(Uniform(0, 1), (T)) 

        for t_ind in range(1, stop = (T-1), step = 1)
            if t_ind == 1
                cut_lb = z_min - mu
                eps_t = trunc_norm_draws(unif_vec[t_ind], 0, sigma, cut_lb)
                zt_vec[t_ind] = mu .+ eps_t
            else
                cut_lb = (z_min - rho * zt_vec[t_ind - 1] - (1 - rho) * mu)
                eps_t = trunc_norm_draws(unif_vec[t_ind], 0, sigma, cut_lb)
                zt_vec[t_ind] = (rho * zt_vec[t_ind - 1] + (1 - rho) * mu .+ eps_t)
            end 
        end 
    end 
    
    default_vec = zeros(Bool, (T))
    c1t_vec = zeros(T)
    c2t_vec = zeros(T)
    Ht_vec = zeros(T)
    wt_vec = zeros(T)
    rt_vec = zeros(T)
    k2t_vec = zeros(T + 1)
    EulErr_vec = zeros(T)
    k2t_vec[1] = k20
    rbart_vec = zeros(T)
    rbart_an_vec = zeros(T)
    t_ind = 1
    default = false

    while (t_ind < (T+1)) && !default
        #print("H_ind = ", H_ind, ",risk_type_ind = ", risk_type_ind, ", risk_val_ind = ", risk_val_ind, ", avgrtp1_ind = ", avgrtp1_ind,
        # ", avgrbart_ind = ", avgrbart_ind, ", S_ind = ", S_ind, ", t_ind = ", t_ind)
        #print(" t_ind = ", t_ind)
        k2t = k2t_vec[t_ind]
        zt = zt_vec[t_ind]
        k2tp1_args = (Hbar, beta, gamma, k20, rho, mu, sigma, x1, nvec,
                      epsilon, alpha, delta, tau, c_min, K_min, A_min,
                      yrs_in_per)
        k2tp1, c1t, Ht, c2t, wt, rt, rbart, rbart_an, default, eul_err = get_k2tp1(k2t, zt, k2tp1_args) 
        k2t_vec[t_ind + 1] = k2tp1[1]
        EulErr_vec[t_ind] = eul_err[1]
        c1t_vec[t_ind] = c1t[1]
        Ht_vec[t_ind] = Ht
        c2t_vec[t_ind] = c2t
        wt_vec[t_ind] = wt
        rt_vec[t_ind] = rt
        rbart_vec[t_ind] = rbart[1]
        rbart_an_vec[t_ind] = rbart_an[1]
        if default
            default_vec[(t_ind+1):end] .= default
        end 
        t_ind = t_ind + 1
    end 

    ut_vec = ((1 - beta) .* log.(c1t_vec[1:(end-1)]) + beta * (1 / (1 - gamma)) .* log.(c2t_vec[2:end].^(1 - gamma)))

    print("H_ind = ", H_ind, ", risk_type_ind = ", risk_type_ind,
          ", risk_val_ind = ", risk_val_ind, ", avgrtp1_ind = ", avgrtp1_ind,
          ", avgrbart_ind = ", avgrbart_ind, ", S_ind = ", S_ind," ")

    return [H_ind, risk_type_ind, risk_val_ind, avgrtp1_ind,
            avgrbart_ind, S_ind, zt_vec, default_vec, c1t_vec,
            c2t_vec, ut_vec, Ht_vec, wt_vec, rt_vec, k2t_vec, rbart_vec,
            rbart_an_vec, EulErr_vec]
end 
