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

"""
Calibrate the parameters
"""
yrs_in_per = 25
beta_an = 0.96
beta = beta_an^yrs_in_per
gamma = 2.2
c_min = 1e-5
K_min = 1e-5

# Firm parameters
alpha = 1/3
epsilon = 1.0  # Inf
delta_an = 0.0
delta = 1 - ((1 - delta_an)^yrs_in_per)
nvec = Array([1.0, 0.0])

# Aggregate shock z parameters
rho_an = 0.95
rho = rho_an^yrs_in_per
mu_an = 0.0
sigma_an = 0.2  # 0.4946
rho_sum = 0.0
rho2_sum = 0.0
for y_ind in range(0,stop = (yrs_in_per-1),step = 1)
    rho_sum = rho_sum + rho_an^y_ind
    rho2_sum = rho2_sum + rho_an^(2*y_ind)
end 
sigma = sqrt(rho2_sum*(sigma_an^2))
mu = mu_an*rho_sum
A_min = 0.0
if A_min == 0.0
    z_min = -Inf
elseif (A_min > 0.0) && (A_min < exp(mu))
    z_min = log(A_min)
elseif A_min >= exp(mu)
    print("Parameter Error: A_min >= e^(mu)")
end

# Set government parameters, transfer parameters, and initial values
Hbar_vec = Array([0.0, 0.05])
# self.Hbar_vec = Array([0.0, 0.05, 0.11, 0.17])
Hbar_size = length(Hbar_vec)
Hbar = Hbar_vec[1]
tau = nothing
z0 = mu

# Set simulation parameters
T = 25
S = 15
rand_seed = 25

# print("sigma")
# print(sigma)
# print(" ")
# print("beta")
# print(beta)
# print(" ")


"""
------------------------------------------------------------------------
Calibrate beta, mu, gamma
------------------------------------------------------------------------
"""

avg_rtp1_size = 3
avg_rtp1_an_vec = Array(range(0.00,stop = 0.04,length = avg_rtp1_size))
avg_Rtp1_vec = (1 .+ avg_rtp1_an_vec).^yrs_in_per
avg_rbart_size = 3
avg_rbart_an_vec = Array(range(-0.02,stop = 0.01, length = avg_rbart_size))
avg_Rbart_vec = (1 .+ avg_rbart_an_vec).^yrs_in_per

# print("avg_Rtp1_vec")
# print(avg_Rtp1_vec)
# print(" ")
# print("avg_Rbart_vec")
# print(avg_Rbart_vec)
# print(" ")

avgRtp1_mat = repeat(reshape(avg_Rtp1_vec,(avg_rtp1_size, 1)), 1, avg_rbart_size)
avgRbart_mat = repeat(reshape(avg_Rbart_vec,(1, avg_rbart_size)), avg_rtp1_size, 1)
avgRtp1_gt_avgRbart = ((avgRtp1_mat - avgRbart_mat) .> 0) 
# print(avgRtp1_gt_avgRbart)

# Calibrate mu using linear production expected MPK
mu_vec = 1.0786 .* fill!(similar(avg_Rtp1_vec), 1)
# mu_vec = log.(avg_Rtp1_vec) .- log(alpha) .- ((sigma^2) / 2)
mu_mat = repeat(reshape(mu_vec,(avg_rtp1_size, 1)), 1, avg_rbart_size)
mu_mat[.!avgRtp1_gt_avgRbart] .= NaN 

# Calibrate beta using Cobb-Douglas expected values expression for beta
beta_vec = (alpha / (1 - alpha)) .* (1 ./ (2 * avg_Rtp1_vec)) 
beta_an_vec = beta_vec.^(1 / yrs_in_per)
beta_mat = repeat(reshape(beta_vec,(avg_rtp1_size, 1)), 1, avg_rbart_size)
beta_mat[.!avgRtp1_gt_avgRbart] .=  NaN

# Calibrate gamma
gamma_mat = ((log.(avgRtp1_mat) - log.(avgRbart_mat)) ./ (sigma^2))
gamma_mat[.!avgRtp1_gt_avgRbart]  .=  NaN

# Calibrate x_1
x1_mat = (((1 - alpha) * exp.(mu_mat .+ ((sigma ^ 2) / 2)) .* ((2 * beta_mat) .^alpha)) .^ (1 / (1 - alpha)))

# Calibrate kbar_2
kbar2_mat = 2 .* beta_mat .* x1_mat
Hbar_vec[2] = 0.05 * kbar2_mat[1, 1] # 1.0786
# print("Hbar_vec = ", Hbar_vec)

# Calibrate x_1 array for different values of x1, given calibration
x1_mat2 = transpose(x1_mat)
x1_arr = zeros(3, avg_rtp1_size, avg_rbart_size) 
x1_arr[:, 1, :] = x1_mat2 
x1_arr[:, 2, :] = 0.5 .* x1_mat2
x1_arr[:, 3, :] = 0.0 .* x1_mat2 

# Calibrate sigma vector of 5% and 10% increases
sigma_vec = zeros(3)
sigma_vec[1] = sigma
sigma_vec[2] = 1.05 * sigma
sigma_vec[3] = 1.10 * sigma

# Calibrate corresponding mu_arr that holds the expected value of the
# TFP shock while expanding the variance. If ExpA is the expected value of
# the TFP shock, then ExpA = exp.(mu .+ (sig^2) / 2), then
# log.(ExpA) = mu .+ (sig^2) / 2
ExpA = mu_mat .+ (sigma ^ 2) / 2
mu_mat2 = transpose(mu_mat)
mu_arr = zeros(3, avg_rtp1_size, avg_rbart_size)
mu_arr[:, 1, :] = mu_mat2
mu_arr[:, 2, :] = mu_mat2
mu_arr[:, 3, :] = mu_mat2
# mu_arr[:, :, 2] = ExpA .- (sigma_vec[2] ^ 2) / 2
# mu_arr[:, :, 3] = ExpA .- (sigma_vec[3] ^ 2) / 2

# print("mu_mat")
# print(mu_mat)
# print(" ")
# print("beta_mat")
# print(beta_mat)
# print(" ")
# print("gamma_mat")
# print(gamma_mat)
# print(" ")
# print("x1_mat")
# print(x1_mat)
# print(" ")
# print("kbar2_mat")
# print(kbar2_mat)
# print(" ")
# print("x1_arr 1")
# print(x1_arr[:, :, 1])
# print(" ")
# print("x1_arr 2")
# print(x1_arr[:, :, 2])
# print(" ")
# print("x1_arr 3")
# print(x1_arr[:, :, 3])
# print(" ")
# print("sigma_vec")
# print(sigma_vec)
# print(" ")
# print("mu_arr 1")
# print(mu_arr[:, :, 1])
# print(" ")
# print("mu_arr 2")
# print(mu_arr[:, :, 2])
# print(" ")
# print("mu_arr 3")
# print(mu_arr[:, :, 3])
# print(" ")

default_arr = zeros(S, T, avg_rbart_size, avg_rtp1_size, 3, 2, Hbar_size)

Random.seed!(rand_seed) 
unif_mat = rand(Uniform(0, 1), (S, T))

# First three dimensions of zt_arr correspond to mu_arr in different order
zt_arr = zeros(S, T, avg_rbart_size, avg_rtp1_size, 3)

cut_lb = 0
eps_t = 0
z_t = 0
z_tm1 = 0

for sig_ind in range(1, stop=(3), step=1)
    sigma = sigma_vec[sig_ind]
    for avgrtp1_ind in range(1, stop=(avg_rtp1_size), step=1)
        for avgrbart_ind in range(1, stop=(avg_rbart_size), step=1)
            mu = mu_arr[avgrbart_ind, sig_ind, avgrtp1_ind]
            for s_ind in range(1, stop=(S), step=1)
                for t_ind in range(1, stop=(T), step=1)
                    unif = unif_mat[s_ind, t_ind]
                    if t_ind == 1 && avgRtp1_gt_avgRbart[avgrtp1_ind,avgrbart_ind] == true
                        cut_lb = z_min - mu
                        eps_t = trunc_norm_draws(unif, 0, sigma, cut_lb)
                        z_t = mu + eps_t
                    elseif ((t_ind > 1) && avgRtp1_gt_avgRbart[avgrtp1_ind,avgrbart_ind] == true)
                        z_tm1 = zt_arr[s_ind, t_ind - 1, avgrbart_ind, avgrtp1_ind, sig_ind]
                        cut_lb = z_min - rho * z_tm1 - (1 - rho) * mu
                        eps_t = trunc_norm_draws(unif, 0, sigma, cut_lb)
                        z_t = rho * z_tm1 + (1 - rho) * mu .+ eps_t 
                    else
                        z_t = NaN
                    end
                zt_arr[s_ind, t_ind, avgrbart_ind,  avgrtp1_ind, sig_ind] = z_t  
                end 
            end 
        end 
    end 
end 

c1t_arr = zero(default_arr)
c2t_arr = zero(default_arr)
ut_arr = zeros(S, T - 1, avg_rbart_size, avg_rtp1_size, 3, 2, Hbar_size)
Ht_arr = zero(default_arr)
wt_arr = zero(default_arr)
rt_arr = zero(default_arr)
k2t_arr = zero(default_arr)
rbart_arr = zero(default_arr)
rbart_an_arr = zero(default_arr)
EulErr_arr = zero(default_arr)
PathTime_arr = zeros(avg_rbart_size, S, avg_rtp1_size, 3, 2, Hbar_size)
s_ind_arr = zeros(avg_rbart_size, S, avg_rtp1_size, 3, 2, Hbar_size)

for rtp1_ind in range(1, stop=(avg_rtp1_size), step=1) 
    for rbart_ind in range(1, stop=(avg_rbart_size), step=1) 
        kbar2_mat2 = (!iszero).(kbar2_mat)
        k2t_arr[:, 1, rbart_ind, rtp1_ind, :, :, :] .=  kbar2_mat2[rtp1_ind, rbart_ind]
    end
end

end


