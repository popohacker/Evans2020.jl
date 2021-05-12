# Import packages 
using Distributions
using Random
using HDF5
using JLD
#à compléter si besoin

include("B19_funcs_NPC.jl") #pour inclure le fichier fonction et que les fonctions soient identifiées par Julia

# Create output directory
cur_path = string(@__DIR__)
output_fldr = "OUTPUT"
output_dir = joinpath(cur_path, output_fldr)
if !ispath(output_dir)
mkpath(output_dir)
end 

yrs_in_per = 25
beta_an = 0.96
beta = beta_an^yrs_in_per
gamma = 2.2
c_min = 1e-5
K_min = 1e-5

# Firm parameters
alpha = 1/3
epsilon = 1.0  # np.inf
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

# Set up parallel processing
max_cores = Sys.CPU_THREADS
print("Cores available on this machine = ", max_cores, " ")
num_workers = min(max_cores, S)
print("Number of workers = ", num_workers)
#client = Client(processes=False) #pas réussi à transcrire cette ligne en julia pour l'instant
#-p num_workers #sort une erreur sur Julia

# print("sigma")
# print(sigma)
# print(" ")
# print("beta")
# print(beta)
# print(" ")

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
avgRtp1_gt_avgRbart = ((avgRtp1_mat - avgRbart_mat) .> 0) #ATTENTION, je ne suis pas certaine qu'ici ce soit bon, dans Python, cela retourne une matrice de true/false tandis que moi j'ai une matrice de 0/1, mais ça me semble bon au vu des résultats suivants
# print(avgRtp1_gt_avgRbart)

# Calibrate mu using linear production expected MPK
mu_vec = 1.0786 .* fill!(similar(avg_Rtp1_vec), 1)
# mu_vec = log.(avg_Rtp1_vec) .- log(alpha) .- ((sigma^2) / 2)
mu_mat = repeat(reshape(mu_vec,(avg_rtp1_size, 1)), 1, avg_rbart_size)
mu_mat[.!avgRtp1_gt_avgRbart] .= NaN #idem peut-être à vérifier mais ça me semble ok

# Calibrate beta using Cobb-Douglas expected values expression for beta
beta_vec = (alpha / (1 - alpha)) .* (1 ./ (2 * avg_Rtp1_vec)) #ATTENTION, le résultat n'est pas exactement identique même si je crois que c'est une erreur d'affichage, car dans Python il me donne 0.25 au premier élément alors que dans Julia il ne l'arrondit pas
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
x1_mat = transpose(x1_mat)
x1_arr = zeros(3, avg_rtp1_size, avg_rbart_size)
x1_arr[:, 1, :] = x1_mat 
x1_arr[:, 2, :] = 0.5 .* x1_mat
x1_arr[:, 3, :] = 0.0 .* x1_mat 

# Calibrate sigma vector of 5% and 10% increases
sigma_vec = zeros(3)
sigma_vec[1] = sigma
sigma_vec[2] = 1.05 * sigma
sigma_vec[3] = 1.10 * sigma

# Calibrate corresponding mu_arr that holds the expected value of the
# TFP shock while expanding the variance. If ExpA is the expected value of
# the TFP shock, then ExpA = exp(mu .+ (sig^2) / 2), then
# log(ExpA) = mu .+ (sig^2) / 2
ExpA = mu_mat .+ (sigma ^ 2) / 2
mu_mat = transpose(mu_mat)
mu_arr = zeros(avg_rtp1_size, avg_rbart_size, 3)
mu_arr[:, 1, :] = mu_mat
mu_arr[:, 2, :] = mu_mat
mu_arr[:, 3, :] = mu_mat
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

#@time begin #à reprendre 

default_arr = zeros(Bool, (Hbar_size, 2, 3, avg_rtp1_size, avg_rbart_size, S, T))

Random.seed!(rand_seed) #pas parfaitement sûre
unif_mat = rand(Uniform(0, 1), (S, T)) #pas parfaitement sûre

# First three dimensions of zt_arr correspond to mu_arr in different order
zt_arr = zeros(3, avg_rtp1_size, avg_rbart_size, S, T)
cut_lb = 0
eps_t = 0
z_t = 0
z_tm1 = 0

for sig_ind in range(1, stop=(3-1), step=1)
    sigma = sigma_vec[sig_ind]
    for avgrtp1_ind in range(1, stop=(avg_rtp1_size-1), step=1)
        for avgrbart_ind in range(1, stop=(avg_rbart_size-1), step=1)
            mu = mu_arr[avgrtp1_ind, avgrbart_ind, sig_ind]
            for s_ind in range(1, stop=(S-1), step=1)
                for t_ind in range(1, stop=(T-1), step=1)
                    unif = unif_mat[s_ind, t_ind]
                    if t_ind == 1 && avgRtp1_gt_avgRbart[avgrtp1_ind,avgrbart_ind]
                        cut_lb = z_min - mu
                        eps_t = trunc_norm_draws(unif, 0, sigma, cut_lb) #plus besoin de mettre funcs devant
                        z_t = mu + eps_t
                    elseif (t_ind > 1) && avgRtp1_gt_avgRbart[avgrtp1_ind,avgrbart_ind]
                        z_tm1 = zt_arr[sig_ind, avgrtp1_ind, avgrbart_ind, s_ind, t_ind - 1] ## REGARDER SI ON BOUGE T_IND - 1 EN PREMIERE PLACE SI ÇA FONCTIONNE
                        cut_lb = z_min - rho * z_tm1 - (1 - rho) * mu
                        eps_t = trunc_norm_draws(unif, 0, sigma, cut_lb)
                        z_t = rho * z_tm1 + (1 - rho) * mu .+ eps_t 
                    else
                        z_t = NaN
                    end
                zt_arr[sig_ind, avgrtp1_ind, avgrbart_ind, s_ind, t_ind] = z_t  
                end 
            end 
        end 
    end 
end 

##### ON CHECK JUSQU'ICI
##### VÉRIFIER L'INDEXATION DANS LA LOOP 
                                                             
c1t_arr = zero(default_arr)
c2t_arr = zero(default_arr)
ut_arr = zeros(Hbar_size, 2, 3, avg_rtp1_size, avg_rbart_size, S, T - 1)
Ht_arr = zero(default_arr)
wt_arr = zero(default_arr)
rt_arr = zero(default_arr)
k2t_arr = zero(default_arr)
rbart_arr = zero(default_arr)
rbart_an_arr = zero(default_arr)
EulErr_arr = zero(default_arr)
PathTime_arr = zeros(Hbar_size, 2, 3, avg_rtp1_size, avg_rbart_size, S)
s_ind_arr = zeros(Hbar_size, 2, 3, avg_rtp1_size, avg_rbart_size, S)

for rtp1_ind in range(1, stop=(avg_rtp1_size-1), step=1)
    for rbart_ind in range(1, stop=(avg_rbart_size-1), step=1)
        k2t_arr = convert.(Float64, k2t_arr)
        k2t_arr[:, :, :, rtp1_ind, rbart_ind, :, 1] .=  kbar2_mat[rtp1_ind, rbart_ind] #à vérifier au vu de ce qu'avait dit Circé la dernière fois mais je pense que c'est bon
        k2t_arr = round.(k2t_arr)
        k2t_arr = convert.(Bool, k2t_arr)
    end
end

for H_ind in 1:Hbar_size
    Hbar_in = Hbar_vec[H_ind]
    for risk_type_ind in 1:2 #0=xval, 1=sigval
        for risk_val_ind in 1:3
            for avgrtp1_ind in 1:avg_rtp1_size
                for avgrbart_ind in 1:avg_rbart_size
                    if avgRtp1_gt_avgRbart[avgrtp1_ind, avgrbart_ind]
                        simulations = []
                        beta_in = beta_mat[avgrtp1_ind, avgrbart_ind]
                        gamma_in = gamma_mat[avgrtp1_ind, avgrbart_ind]
                        k20_in = kbar2_mat[avgrtp1_ind, avgrbart_ind]
                        if risk_type_ind == 1
                            mu_in = mu_mat[avgrtp1_ind, avgrbart_ind]
                            sigma_in = sigma_vec[1]
                            x1_in = x1_arr[avgrtp1_ind, avgrbart_ind,
                                           risk_val_ind]
                            z0_vec_in = zt_arr[1, avgrtp1_ind,
                                           avgrbart_ind, :, 1] ##pas sure
                        elseif risk_type_ind == 2 
                            mu_in = mu_arr[avgrtp1_ind, avgrbart_ind,
                                           risk_val_ind]
                            sigma_in = sigma_vec[risk_val_ind]
                            x1_in = x1_mat[avgrtp1_ind, avgrbart_ind]
                            z0_vec_in = zt_arr[risk_val_ind,
                                               avgrtp1_ind,
                                               avgrbart_ind, :, 1]## pas sure
                        end

                        for s_ind in 1:S
                            z0_in = z0_vec_in[s_ind]
                            if risk_type_ind == 1
                                zt_vec_in = zt_arr[1, avgrtp1_ind, avgrbart_ind,
                                           s_ind, :]
                            elseif risk_type_ind == 2
                                zt_vec_in = zt_arr[risk_type_ind, avgrtp1_ind,
                                           avgrbart_ind, s_ind, :]
                            end
                            timepaths_s = sim_timepath(
                                Hbar_in, beta_in, gamma_in, k20_in,
                                sigma_in, x1_in, T, z0_in, z_min, rho,
                                mu_in, nvec, epsilon, alpha, delta, tau,
                                c_min, K_min, A_min, yrs_in_per,
                                H_ind= H_ind,
                                risk_type_ind=risk_type_ind,
                                risk_val_ind=risk_val_ind,
                                avgrtp1_ind=avgrtp1_ind,
                                avgrbart_ind=avgrbart_ind, S_ind=s_ind,
                                zt_vec=zt_vec_in,
                                rand_seed=rand_seed)
                            append!(simulations, timepaths_s)
                        end #not sure if the for loop ends here or later
                        
                        #Original line is : simulations = delayed(simulations).compute()
                        #I did not add this line because if we don't use delayed I think it's not needed 

                        for s_ind in 1:S
                            s_ind_arr[H_ind, risk_type_ind,
                                      risk_val_ind, avgrtp1_ind,
                                      avgrbart_ind, s_ind] = simulations[s_ind][6] # original S_ind
                            default_arr[H_ind, risk_type_ind,
                                      risk_val_ind, avgrtp1_ind,
                                      avgrbart_ind, s_ind, :] = simulations[s_ind][8]  # default_vec
                            c1t_arr[H_ind, risk_type_ind, risk_val_ind,
                                      avgrtp1_ind, avgrbart_ind, s_ind,
                                      :] = simulations[s_ind][9]  # c1t_vec
                            c2t_arr[H_ind, risk_type_ind, risk_val_ind,
                                      avgrtp1_ind, avgrbart_ind, s_ind,
                                      :] = simulations[s_ind][10]  # c2t_vec
                            ut_arr[H_ind, risk_type_ind, risk_val_ind,
                                     avgrtp1_ind, avgrbart_ind, s_ind,
                                     :] = simulations[s_ind][11]  # ut_vec
                            Ht_arr[H_ind, risk_type_ind, risk_val_ind,
                                     avgrtp1_ind, avgrbart_ind, s_ind,
                                     :] = simulations[s_ind][12]  # Ht_vec
                            wt_arr[H_ind, risk_type_ind, risk_val_ind,
                                     avgrtp1_ind, avgrbart_ind, s_ind,
                                     :] = simulations[s_ind][13]  # wt_vec
                            rt_arr[H_ind, risk_type_ind, risk_val_ind,
                                     avgrtp1_ind, avgrbart_ind, s_ind,
                                     :] = simulations[s_ind][14]  # rt_vec
                            k2t_arr[H_ind, risk_type_ind, risk_val_ind,
                                      avgrtp1_ind, avgrbart_ind, s_ind,
                                      :] = pop!(simulations[s_ind][15])  # k2t_vec[:-1]
                            rbart_arr[H_ind, risk_type_ind,
                                        risk_val_ind, avgrtp1_ind,
                                        avgrbart_ind, s_ind, :] = simulations[s_ind][16]  # rbart_vec
                            rbart_an_arr[H_ind, risk_type_ind,
                                           risk_val_ind, avgrtp1_ind,
                                           avgrbart_ind, s_ind, :] = simulations[s_ind][17]  # rbart_an_vec
                            EulErr_arr[H_ind, risk_type_ind,
                                         risk_val_ind, avgrtp1_ind,
                                         avgrbart_ind, s_ind, :] = simulations[s_ind][18]  # EulErr_vec
                            PathTime_arr[H_ind, risk_type_ind,
                                           risk_val_ind, avgrtp1_ind,
                                           avgrbart_ind, s_ind] = simulations[s_ind][19]  # path_time
                        end
                    else  # avg_Rtp1 <= avg_rbart
                        s_ind_arr[H_ind, risk_type_ind, risk_val_ind,
                        avgrtp1_ind, avgrbart_ind, :] .= NaN 
                        default_arr[H_ind, risk_type_ind, risk_val_ind,
                                    avgrtp1_ind, avgrbart_ind, :, :] .= NaN # default_vec
                        c1t_arr[H_ind, risk_type_ind, risk_val_ind,
                                avgrtp1_ind, avgrbart_ind, :, :] .= NaN # c1t_vec
                        c2t_arr[H_ind, risk_type_ind, risk_val_ind,
                                avgrtp1_ind, avgrbart_ind, :, :] .= NaN # c2t_vec
                        ut_arr[H_ind, risk_type_ind, risk_val_ind,
                               avgrtp1_ind, avgrbart_ind, :, :] .= NaN  # ut_vec
                        Ht_arr[H_ind, risk_type_ind, risk_val_ind,
                               avgrtp1_ind, avgrbart_ind, :, :] .= NaN  # Ht_vec
                        wt_arr[H_ind, risk_type_ind, risk_val_ind,
                               avgrtp1_ind, avgrbart_ind, :, :] .= NaN # wt_vec
                        rt_arr[H_ind, risk_type_ind, risk_val_ind,
                               avgrtp1_ind, avgrbart_ind, :, :] .= NaN # rt_vec
                        k2t_arr[H_ind, risk_type_ind, risk_val_ind,
                                avgrtp1_ind, avgrbart_ind, :, :] .= NaN # k2t_vec[:-1]
                        rbart_arr[H_ind, risk_type_ind, risk_val_ind,
                                  avgrtp1_ind, avgrbart_ind, :, :] .= NaN # rbart_vec
                        rbart_an_arr[H_ind, risk_type_ind, risk_val_ind,
                                     avgrtp1_ind, avgrbart_ind, :,
                                     :] .= NaN  # rbart_an_vec
                        EulErr_arr[H_ind, risk_type_ind, risk_val_ind,
                                   avgrtp1_ind, avgrbart_ind, :, :] .= NaN # EulErr_vec
                        PathTime_arr[H_ind, risk_type_ind, risk_val_ind,
                                     avgrtp1_ind, avgrbart_ind, :] .= NaN  # path_time
                    end
                end
            dict_endog_new = Dict(
                "unif_mat" => unif_mat,
                "zt_arr" => zt_arr ,
                "c1t_arr" =>  c1t_arr,
                "c2t_arr"=>  c2t_arr,
                "ut_arr" =>  ut_arr,
                "Ht_arr" =>  Ht_arr,
                "wt_arr" =>  wt_arr,
                "rt_arr" =>  rt_arr,
                "rbart_arr" =>  rbart_arr,
                "rbart_an_arr" => rbart_an_arr,
                "k2t_arr" =>  k2t_arr,
                "EulErr_arr" =>  EulErr_arr,
                "PathTime_arr" =>  PathTime_arr,
                "default_arr" =>  default_arr,
                "s_ind_arr" =>  s_ind_arr)
            
            #exec('outputfile = os.path.join(output_dir, \'dict_endog_' +
                #str(H_ind) + str(risk_type_ind) + str(risk_val_ind) +
                #'.pkl\')')
            #exec('pickle.dump(dict_endog_new, open(outputfile, ' +
                #'\'wb\'))')  
                
            save("dict_endog_$(H_ind)$(risk_type_ind)$(risk_val_ind).jld", "dict_endog_new", dict_endog_new)   
            #ça save dans le directory dans lequel on est -- s'assurer qu'on est dans output
            end
        end
    end
end
    
default_p1 = cat(6,zeros(Bool, (H_ind, risk_type_ind, risk_val_ind, avgrtp1_ind,
                        avgrbart_ind, S, 1)),
              default_arr[:, :, :, :, :, :, 1]) #pas fini

zt_arr = zeros(3, avg_rtp1_size, avg_rbart_size, S, T)
size(zt_arr)
zt_arr_macro = repeat(reshape(zt_arr,(1, 1, 3, avg_rtp1_size,
                                       avg_rbart_size, S, T)),
                       Hbar_size, 2, 1, 1, 1, 1, 1)

Kt_arr = (1 - default_p1) .* k2t_arr
Y_args = (nvec, epsilon, alpha)
Yt_arr = (1 - default_p1) .* get_Y(Kt_arr, zt_arr_macro, Y_args)
Ct_arr = (1 - default_p1) .* get_C(c1t_arr, c2t_arr)

dict_params = Dict(
        "yrs_in_per" => yrs_in_per,
        "beta_an" => beta_an,
        "beta" => beta,
        "gamma" => gamma,
        "c_min" => c_min,
        "K_min" => K_min,
        "nvec" => nvec,
        "n1" => nvec[0],
        "n2" => nvec[1],
        "alpha" => alpha,
        "epsilon" => epsilon,
        "delta_an" => delta_an,
        "delta" => delta,
        "rho_an" => rho_an,
        "rho" => rho,
        "mu_an" => mu_an,
        "sigma_an" => sigma_an,
        "sigma" => sigma,
        "mu" => mu,
        "A_min" => A_min,
        "z_min" => z_min,
        "Hbar_vec" => Hbar_vec,
        "Hbar_size" => Hbar_size,
        "Hbar" => Hbar,
        "tau" => tau,
        "T" => T,
        "S" => S,
        "rand_seed" => rand_seed,
        "max_cores" => max_cores,
        "num_workers" => num_workers,
        "avg_rtp1_size" => avg_rtp1_size,
        "avg_rtp1_an_vec" => avg_rtp1_an_vec,
        "avg_Rtp1_vec" => avg_Rtp1_vec,
        "avg_rbart_size" => avg_rbart_size,
        "avg_rbart_an_vec" => avg_rbart_an_vec,
        "avg_Rbart_vec" => avg_Rbart_vec,
        "avgRtp1_mat" => avgRtp1_mat,
        "avgRbart_mat" => avgRbart_mat,
        "avgRtp1_gt_avgRbart" => avgRtp1_gt_avgRbart,
        "mu_vec" => mu_vec,
        "mu_mat" => mu_mat,
        "mu_arr" => mu_arr,
        "beta_vec" => beta_vec,
        "beta_mat" => beta_mat,
        "gamma_mat" => gamma_mat,
        "x1_mat" => x1_mat,
        "x1_arr" => x1_arr,
        "kbar2_mat" => kbar2_mat,
        "sigma_vec" => sigma_vec,
        "ExpA" => ExpA)

    dict_endog = Dict(
        "unif_mat" => unif_mat,
                "zt_arr" => zt_arr ,
                "c1t_arr" =>  c1t_arr,
                "c2t_arr"=>  c2t_arr,
                "ut_arr" =>  ut_arr,
                "Ht_arr" =>  Ht_arr,
                "wt_arr" =>  wt_arr,
                "rt_arr" =>  rt_arr,
                "rbart_arr" =>  rbart_arr,
                "rbart_an_arr" => rbart_an_arr,
                "k2t_arr" =>  k2t_arr,
                "EulErr_arr" =>  EulErr_arr,
                "PathTime_arr" =>  PathTime_arr,
                "Kt_arr" => Kt_arr,
                "Yt_arr"=> Yt_arr,
                 "Ct_arr" => Ct_arr,
                "default_arr" =>  default_arr,
                "s_ind_arr" =>  s_ind_arr)
                #"total_time" => total_time) #à reprendre

                #ca
results_sims = Dict("dict_params" => dict_params, "dict_endog" => dict_endog)
#outputfile = os.path.join(output_dir, "results_sims.pkl") #à reprendre

#ca g as fe #à reprendre
#pickle.dump(results_sims, open(outputfile, "wb")) #lol wat 
#pickle.dump(pythonObject, pickleDestination, pickle_protocol=None, *, fix_imports=True) #j'ai pas cette ligne moi sur mon fichier python

save("results_sims.jld", "results_sims", results_sims)
    
#end #celui de time
