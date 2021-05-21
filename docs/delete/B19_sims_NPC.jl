#Import packages 
#à compléter

#Create output directory
#à compléter


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
x1_arr = zeros(avg_rtp1_size, avg_rbart_size, 3)
x1_arr[:, :, 1] = x1_mat
x1_arr[:, :, 2] = 0.5 * x1_mat
x1_arr[:, :, 3] = 0.0 * x1_mat #ATTENTION, petit doute sur si le résultat correspond, la présentation des matrices est différentes

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
mu_arr = zeros(avg_rtp1_size, avg_rbart_size, 3)
mu_arr[:, :, 1] = mu_mat
mu_arr[:, :, 2] = mu_mat
mu_arr[:, :, 3] = mu_mat
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
