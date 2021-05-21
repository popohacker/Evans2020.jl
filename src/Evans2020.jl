module Evans2020
using JLD2
using FileIO
using Statistics
using Distributions
using Random
using LinearAlgebra
using Roots
using SymEngine
using QuadGK
using DifferentialEquations
using CommonSolve
using Optim
using LeastSquaresOptim 


"""
trunc_norm_draws(unif_vals::Any, mu::Any, sigma::Any, cut_lb::Any = nothing, cut_ub::Any = nothing)
    if cut_lb == nothing && cut_ub == nothing

    --------------------------------------------------------------------
    Draw (N x S) matrix of random draws from a truncated normal
    distribution based on a normal distribution with mean mu and
    standard deviation sigma and cutoffs (cut_lb, cut_ub). These draws
    correspond to an (N x S) matrix of randomly generated draws from a
    uniform distribution U(0,1).
    --------------------------------------------------------------------
    INPUTS:
    unif_vals = (N, S) matrix, (N,) vector, or scalar in (0,1), random
                draws from uniform U(0,1) distribution
    mu        = scalar, mean of the nontruncated normal distribution
                from which the truncated normal is derived
    sigma     = scalar > 0, standard deviation of the nontruncated
                normal distribution from which the truncated normal is
                derived
    cut_lb    = scalar or string, ='None' if no lower bound cutoff is
                given, otherwise is scalar lower bound value of
                distribution. Values below this cutoff have zero
                probability
    cut_ub    = scalar or string, ='None' if no upper bound cutoff is
                given, otherwise is scalar lower bound value of
                distribution. Values below this cutoff have zero
                probability

    OTHER FUNCTIONS AND FILES CALLED BY THIS FUNCTION:
        scipy.stats.norm()

    OBJECTS CREATED WITHIN FUNCTION:
    cut_ub_cdf  = scalar in [0, 1], cdf of N(mu, sigma) at upper bound
                  cutoff of truncated normal distribution
    cut_lb_cdf  = scalar in [0, 1], cdf of N(mu, sigma) at lower bound
                  cutoff of truncated normal distribution
    unif2_vals  = (N, S) matrix, (N,) vector, or scalar in (0,1),
                  rescaled uniform derived from original.
    tnorm_draws = (N, S) matrix, (N,) vector, or scalar in (0,1),
                  values drawn from truncated normal PDF with base
                  normal distribution N(mu, sigma) and cutoffs
                  (cut_lb, cut_ub)

    FILES CREATED BY THIS FUNCTION: None

    RETURNS: tnorm_draws
    --------------------------------------------------------------------


"""
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


"""
    get_Y(k2t, zt, args)

    --------------------------------------------------------------------
    Calculate aggregate output
    --------------------------------------------------------------------
    INPUTS:
    Kt   =
    Lt   =
    zt   =
    args =

    RETURNS: Yt
"""

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

"""
sim_timepath(
    Hbar, beta, gamma, k20, sigma, x1, T, z0, z_min, rho, mu, nvec,
    epsilon, alpha, delta, tau, c_min, K_min, A_min, yrs_in_per,
    H_ind, risk_type_ind, risk_val_ind,
    avgrtp1_ind, avgrbart_ind, S_ind, zt_vec,
    rand_seed)

Calls the arguments to run the simulation

"""

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



export      sim_timepath,
            trunc_norm_draws,
            get_Y

end


