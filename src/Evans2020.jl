module Evans2020
using JLD2
using FileIO
using Statistics
using Distributions
using Random
using HDF5
using LinearAlgebra
using Roots
using SymEngine
using QuadGK
using DifferentialEquations
using CommonSolve
using Optim
using LeastSquaresOptim 



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






end


