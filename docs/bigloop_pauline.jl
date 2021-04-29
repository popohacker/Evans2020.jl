using HDF5
using JLD 


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


