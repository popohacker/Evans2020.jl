using JLD2
using FileIO
using Statistics

#change directory, go to our output directory
cd("C:\\Users\\maillet_c\\117321-V1\\code\\e1_mucnst")

# Load pickled results
cur_dir = string(@__DIR__)
output_dir = joinpath(cur_dir, "OUTPUT") 

#Je genere les arrays qu'on devrait avoir genere a la fin de bigloop et je les pickle
#for reference 
#ut_arr = zeros(15, 25, 3, 3, 3, 2, 2)
#rbart_an_arr = zero(default_arr)
#default_arr = zeros(Bool, (15, 25, 3, 3, 3, 2, 2))

for H_ind in range(1, stop=2, step=1) 
    for risk_type_ind in range(1, stop=2, step=1) 
        for risk_val_ind in range(1, stop=3, step=1)
            ut_arr = ones(15, 24, 3, 3, risk_val_ind , risk_type_ind, H_ind) #I create an array of one so that in 2nd loop I can divide by its values
            rbart_an_arr = ones(15, 25, 3, 3, risk_val_ind , risk_type_ind, H_ind)
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
            global ut_arr = load("dict_endog_$h$r1$r2.jld2", "ut_arr") 
            @eval $(Symbol("ut_arr_$h$r1$r2")) = ut_arr[ :, :, :, :, r2, r1, h]
            print("ut_arr_$h$r1$r2")
            global rbart_an_arr = load("dict_endog_$h$r1$r2.jld2", "rbart_an_arr") 
            @eval $(Symbol("rbart_an_arr_$h$r1$r2")) = rbart_an_arr[ :, :, :, :, r2, r1, h]                   
        end
     end
end

# Solve for percent difference in average welfare matrices
avg_rtp1_size = 3
avg_rbart_size = 3

for risk_type_ind in range(1, stop=2, step=1)
  for risk_val_ind in range(1, stop=3, step=1)
    # ut_pctdif_1?? = np.zeros((?, ?))
    global r1  =  risk_type_ind  
    global r2 = risk_val_ind
    @eval $(Symbol("ut_pctdif_1$r1$r2")) = zeros((avg_rtp1_size, avg_rbart_size))
      for avgrtp1_ind in range(1,stop = (avg_rtp1_size), step = 1) 
        for avgrbart_ind in range(1,stop = (avg_rbart_size), step = 1)
          global t = avgrtp1_ind
          global b  =  avgrbart_ind 
          eval(Meta.parse("ut_mat_1$r1$r2 = ut_arr_1$r1$r2[ :, :, $b, $t]"))
          eval(Meta.parse("ut_mat_2$r1$r2 = ut_arr_2$r1$r2[ :, :, $b, $t]"))
          eval(Meta.parse("avg_ut_1$r1$r2 = mean(ut_mat_1$r1$r2[.!isnan.(ut_mat_1$r1$r2)])")) 
          eval(Meta.parse("avg_ut_2$r1$r2 = mean(ut_mat_2$r1$r2[.!isnan.(ut_mat_2$r1$r2)])")) 
          eval(Meta.parse("ut_pctdif_1$r1$r2[b, t] = (avg_ut_2$r1$r2 - avg_ut_1$r1$r2)/avg_ut_1$r1$r2")) 
          println("ut_pctdif_1$r1$r2 for Cobb-Douglas, mu variable")
          eval(Meta.parse("println(ut_pctdif_1$r1$r2)")) #amelioration possible: ne pas applatir les rows (rows sont displayes les uns a cote des autres)
          println()
        end
    end
  end
end
