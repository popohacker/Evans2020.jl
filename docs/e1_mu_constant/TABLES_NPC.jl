using JLD2
using FileIO
using Statistics

#set directory to output, where B19_sims generated pickled output
cd("-")
cur_dir = string(@__DIR__)
output_dir = joinpath(cur_dir, "OUTPUT") 


# Load pickled results
for H_ind in range(1, stop=2, step=1) 
    for risk_type_ind in range(1, stop=2, step=1) 
        for risk_val_ind in range(1, stop=3, step=1)     
            global h = H_ind
            global r1  =  risk_type_ind
            global r2 = risk_val_ind
            println("dict_endog_$h$r1$r2")
            global ut_arr = load("dict_endog_$h$r1$r2.jld2", "ut_arr") 
            eval(Meta.parse("ut_arr_$h$r1$r2 = ut_arr[ :, :, :, :, r2, r1, h]")) 
            println("ut_arr_$h$r1$r2")
            global rbart_an_arr = load("dict_endog_$h$r1$r2.jld2", "rbart_an_arr") 
            eval(Meta.parse("rbart_an_arr_$h$r1$r2 = rbart_an_arr[ :, :, :, :, r2, r1, h]"))                  
         end
     end
 end

# Solve for percent difference in average welfare matrices
avg_rtp1_size = 3
avg_rbart_size = 3

#Generate Paper Table 2
for risk_type_ind in range(1, stop=2, step=1)
  for risk_val_ind in range(1, stop=3, step=1)
     global r1 = risk_type_ind  
     global r2 = risk_val_ind
     eval(Meta.parse("ut_pctdif_1$r1$r2 = zeros((avg_rbart_size, avg_rtp1_size))"))
     for avgrtp1_ind in 1:3 
       for avgrbart_ind in 1:3
         global t = avgrtp1_ind
         global b = avgrbart_ind 
         eval(Meta.parse("ut_mat_1$r1$r2 = ut_arr_1$r1$r2[ :, :, $b, $t]"))
         eval(Meta.parse("ut_mat_2$r1$r2 = ut_arr_2$r1$r2[ :, :, $b, $t]"))
         eval(Meta.parse("avg_ut_1$r1$r2 = mean(ut_mat_1$r1$r2[.!isnan.(ut_mat_1$r1$r2)])")) 
         eval(Meta.parse("avg_ut_2$r1$r2 = mean(ut_mat_2$r1$r2[.!isnan.(ut_mat_2$r1$r2)])")) 
         eval(Meta.parse("ut_pctdif_1$r1$r2[b, t] = (avg_ut_2$r1$r2 - avg_ut_1$r1$r2)/avg_ut_1$r1$r2")) 
        end
      end
      println("ut_pctdif_1$r1$r2 for Cobb-Douglas, mu variable")
      eval(Meta.parse("println(transpose(ut_pctdif_1$r1$r2))")) 
      println()
   end
end

