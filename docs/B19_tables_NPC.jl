# Import packages
using JLD2
using HDF5
using FileIO

#change directory, go to our output directory
cd("-")

# Load pickled results
cur_dir = string(@__DIR__)
output_dir = joinpath(cur_dir, "OUTPUT") 

#Je genere les arrays qu'on devrait aoir genere a la fin de bigloop et je les pickle
#for reference 
#ut_arr = zeros(15, 25, 3, 3, 3, 2, 2)
#rbart_an_arr = zero(default_arr)
#default_arr = zeros(Bool, (15, 25, 3, 3, 3, 2, 2))

for H_ind in range(1, stop=2, step=1) 
    for risk_type_ind in range(1, stop=2, step=1) 
        for risk_val_ind in range(1, stop=3, step=1)
            ut_arr = zeros(15, 25, 3, 3, risk_val_ind , risk_type_ind, H_ind)
            default_arr = zeros(Bool, (15, 25, 3, 3, risk_val_ind , risk_type_ind, H_ind))
            rbart_an_arr = zero(default_arr)
         save("dict_endog_$(H_ind)$(risk_type_ind)$(risk_val_ind).jld2", Dict("ut_arr" =>  ut_arr,"rbart_an_arr" => rbart_an_arr))
        end
     end
end

#boom je les unpickle, jleur donne meme des noms etou
#ici le code de evans reprends

for H_ind in range(1, stop=2, step=1) 
    for risk_type_ind in range(1, stop=2, step=1) 
        for risk_val_ind in range(1, stop=3, step=1)     
            global h = H_ind
            global r1  =  risk_type_ind  
            global r2 = risk_val_ind
            filename = load("dict_endog_$h$r1$r2.jld2")
            global ut_arr = load("dict_endog_$h$r1$r2.jld2", "ut_arr") #ca c bon
            #run("ut_arr_$(H_ind)$(risk_type_ind)$(risk_val_ind)" = ut_arr[$(H_ind),$(risk_type_ind),$(risk_val_ind), :, :, :, :]) #ca ne marchera pas, a modifier, notamment le nom de ce que j'appelle
            @eval $(Symbol("ut_arr_$h$r1$r2")) = ut_arr[ :, :, :, :, r2, r1, h]
            print("ut_arr_$h$r1$r2")
            global rbart_an_arr = load("dict_endog_$h$r1$r2.jld2", "rbart_an_arr") #ca ca devrait aller
            @eval $(Symbol("rbart_an_arr_$h$r1$r2")) = rbart_an_arr[ :, :, :, :, r2, r1, h] #ca ne marchera pas, a modifier                    
        end
     end
end

#je genere des dict pour les acceder dans la boucle d'apres
ut_arr_1 = Dict(:ut_arr_111 => ut_arr_111, :ut_arr_112 => ut_arr_112, :ut_arr_113 => ut_arr_113, :ut_arr_121 => ut_arr_121, :ut_arr_122 => ut_arr_122, :ut_arr_123 => ut_arr_123)
ut_arr_2 = Dict(:ut_arr_211 => ut_arr_211, :ut_arr_212 => ut_arr_212, :ut_arr_213 => ut_arr_213, :ut_arr_221 => ut_arr_221, :ut_arr_222 => ut_arr_222, :ut_arr_223 => ut_arr_223)
#ut_arr_2[:ut_arr_211][1, :, :, :, :, :, :] l'idee c'est de remplacer 11 par nos strings

# Solve for percent difference in average welfare matrices
avg_rtp1_size = 3
avg_rbart_size = 3

for risk_type_ind in range(1, stop=2, step=1)
    for risk_val_ind in range(1, stop=3, step=1)
        # ut_pctdif_1?? = np.zeros((?, ?))
        global r1  =  risk_type_ind  
        global r2 = risk_val_ind
        @eval $(Symbol("ut_pctdif_1$r1$r2") = zeros((avg_rtp1_size, avg_rbart_size))
        for avgrtp1_ind in range(1,stop = (avg_rtp1_size), step = 1) 
          for avgrbart_ind in range(1,stop = (avg_rbart_size), step = 1)
            global t = avgrtp1_ind
            global b  =  avgrbart_ind  
              # ut_mat_0?? = ut_arr_0??[?, ?, :, :]
              #@eval $(Symbol("ut_mat_1$r1$r2") = ut_arr_1$r1$r2[ :, :, b, t]) #probablement un pb 
              #@eval $(Symbol("ut_mat_2$r1$r2") = "ut_arr_2$r1$r2"[ :, :, b, t]) #probablement un pb 
              #global ut_arr_1 = ut_arr_1$r1$r2
              @eval $(Symbol("ut_mat_1$r1$r2") = ut_arr_1[ :, :, b, t]) #probablement un pb ut_arr_2[:ut_arr_211][1, :, :, :, :, :, :]
              #@eval $(Symbol("ut_mat_2$r1$r2") = ut_arr_2[ :, :, b, t]) #probablement un pb 
              #@eval $(Symbol("avg_ut_1$r1$r2") = mean("ut_mat_1$r1$r2"[.!isnan(ut_mat_1$r1$r2])) #probablement un pb    
              #@eval $(Symbol("avg_ut_2$r1$r2") = mean("ut_mat_2$r1$r2"[.!isnan(ut_mat_2$r1$r2])) #probablement un pb 
              #@eval $(Symbol("ut_pctdif_1$r1$r2[b, t]") = "avg_ut_2$r1$r2- avg_ut_1$r1$r2/avg_ut_1$r1$r2"))) #probablement un pb 
              #print("ut_pctdif_1$r1$r2 for Cobb-Douglas, mu variable") #c bon
              #print("ut_pctdif_1$r1$r2") #c bon

              #global rbart_an_arr = load("dict_endog_$h$r1$r2.jld2", "rbart_an_arr") #ca ca devrait aller
              #@eval $(Symbol("rbart_an_arr_$h$r1$r2")) = rbart_an_arr[ :, :, :, :, r2, r1, h] #ca ne marchera pas, a modifier  
          end
        end
     end
end 

#plus petite boucle pour tester les command une a une
ut_arr_1 = dict()
for risk_type_ind in range(1, stop=2, step=1)
    for risk_val_ind in range(1, stop=3, step=1)
        global r1  =  risk_type_ind  
        global r2 = risk_val_ind
        @eval $(Symbol("ut_pctdif_1$r1$r2")) = zeros((avg_rtp1_size, avg_rbart_size))
        for avgrtp1_ind in range(1,stop = 4, step = 1) 
          for avgrbart_ind in range(1,stop = 4, step = 1)
            global t = avgrtp1_ind
            global b  =  avgrbart_ind  
            global ut_arr_1 = ut_arr_1$r1$r2
            @eval $(Symbol("ut_mat_1$r1$r2")) = ut_arr_1[:ut_arr_2$r1$r2][ :, :, b, t] #probablement un pb 
          end
        end
     end
end 

