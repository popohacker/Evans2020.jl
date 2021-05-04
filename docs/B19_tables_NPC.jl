# Import packages
#à compléter si nécessaire

# Load pickled results
cur_dir = string(@__DIR__)
output_dir = joinpath(cur_dir, "OUTPUT")

for H_ind in range(0, stop=1, step=1) #attention, peut-être qu'il faudra changer le range pour qu'il commence par 1 --> à voir par la suite
    for risk_type_ind in range(0, stop=1, step=1) #idem
        for risk_val_ind in range(0, stop=2, step=1) #idem      
            filename = string("dict_endog_$(H_ind)$(risk_type_ind)$(risk_val_ind)")
            print("filename= ", filename)
            #run(string(filename, "_path = joinpath(output_dir, '", filename , ".pkl')")) pas terminé 
      
      
      end 
    end 
end 
