##Kelly Clark
##08/04/2021
##permutation testing: inter/intra coil volumetric differences Siemens Sites 

library(tidyverse)
library(reshape2)
library(dplyr)

 #BWH
 subjects_1 = c("Subject","list","for","BWH")
 #Penn
 subjects_2 = c("Subject","list","for","Penn")

df = read.csv("./path/to/df/containing/volumetric info")
df = melt(df)

#filter out NIH since recieve coil wasn't an issue here
BWH_Penn_df = df %>% filter(Site!="NIH")
#Distortion Corrected
BWH_Penn_D = BWH_Penn_df %>% filter(Type=="D")
#non-Distortion Corrected
BWH_Penn_ND = BWH_Penn_df %>% filter(Type!="D")

#Mean or Median (change file path accordingly) inter/intra ratios of observed data
D_actual_ratios = read.csv("./path/to/distortion-corrected/ratios")
ND_actual_ratios = read.csv("./path/to/non-distortion-corrected/ratios.csv")

#exclude Atropos WM bc outliers
struc_list = c("JLF.Thal.Vol","First.Thal.Vol","Mimosa.Lesion.Vol",
               "JLF.White.Matter.Vol","JLF.Grey.Matter.Vol","FAST.TBV","Atropos.Grey.Matter.Vol")
site_list = c("BWH","Penn")
set.seed(1979)

#initialize empty df to store pvals
permtest_inter_intra_coil_type_Distortion_Corrected_pvals_means = data.frame()

#this is a chaotic way to do this (bunch of for loops) but it works lol!
#loop through structures 
   for(k in 1:length(struc_list)){
      Structure_1 = struc_list[k]  
        #subset data frame by structure[k]
        struc_df = BWH_Penn_D %>% filter(variable==Structure_1)
      #get mean (or median) ratios from observed data -->used later on to calculate pvals
      inter_intra_ratio_actual_data_df = D_actual_ratios %>% filter(Structure==Structure_1)
      inter_intra_actual_ratio = inter_intra_ratio_actual_data_df$ratio_inter_intra_avg_coil_diffs
      
      #number of permutations
      P = 10000
      
      #empty list to store some things...not necessary, but I used these to check things
      list_of_ratios = vector()
      list_of_ratio_comparisons = vector()
      list_of_perm_fails = vector()
      
      #loop through all 10000 permutations...yes,really. Very chaotic and bad way to do this :)
      for(l in 1:P){
        #more empty lists to store things :')
        coil_inter_list = vector()
        coil_intra_list = vector()
       
        #surprise! another for loop!  
        for (j in 1:length(site_list)){
          Site1 = site_list[j]
          vols_df = struc_df %>% filter(Site==Site1)
          n = length(vols_df$Site)
          variable = vols_df$Coil
          
          PermSamples = matrix(0, nrow=n, ncol=P)
          #shuffle coil labels
          PermSamples[,l] = sample(variable, size=n, replace=FALSE)
          PermDF = as.data.frame(PermSamples[,l])
          
          #add permuted coil labels to df containing volumetric info
          vols_df$Coil = PermDF$`PermSamples[, l]`
         
          #if site=Penn, use subjects_2 since one subjcet was not imaged at Penn
           if (Site1=="Penn"){
            # print("ok")
          for (i in 1:length(subjects_2)){
            Subject = subjects_2[i]
            
            #if coil type is same within pair, calculate abs diff and append to intra coil diffs list
            pair = vols_df %>% filter(Subj==Subject)
            if (pair[,4,1][1]==pair[,4,1][2]){
              intracoil_diffs = dist(pair$value)
              coil_intra_list = append(coil_intra_list, intracoil_diffs)
              
            }else{
              #if coil type is different within pair, calculate abs diff and append to inter coil diffs list
              
              intercoil_diffs = dist(pair$value)

              coil_inter_list = append(coil_inter_list, intercoil_diffs)
              
               }

          }
          } else{ for (i in 1:length(subjects_1)){
            Subject = subjects_1[i]
            
            pair = vols_df %>% filter(Subj==Subject)
            if (pair[,4,1][1]==pair[,4,1][2]){
              intracoil_diffs = dist(pair$value)
              coil_intra_list = append(coil_intra_list, intracoil_diffs)
              
            }else{
              intercoil_diffs = dist(pair$value)
              coil_inter_list = append(coil_inter_list, intercoil_diffs)
              
            }
            
          }
          }
        }
        #if permuted data contains no intercoil pairs, move onto next permutation and keep tally of # of times this happens
        message_1 = paste0("Permutation fail: ", Structure_1)
        if (length(coil_inter_list)==0){
          x="fail"
         list_of_perm_fails = append(list_of_perm_fails,x)
          cat(message_1)
         next  
        #else, calculate mean (or median) or inter and intra coil diffs list
          }else{
        avg_coil_intra_both_sites = mean(coil_intra_list)
        avg_coil_inter_both_sites = mean(coil_inter_list)
        
        #calculate inter/intra ratio
        ratio_inter_intra_avg_coil_diffs = avg_coil_inter_both_sites/avg_coil_intra_both_sites
        
        #add to list of ratios...this is pointless honestly, I initially used it to check that my code was working
        list_of_ratios = append(list_of_ratios, ratio_inter_intra_avg_coil_diffs)
        
        #if ratio from permuted data >= ratio from observed data, make note of this 
        if (ratio_inter_intra_avg_coil_diffs>=inter_intra_actual_ratio){ 
          Ratio_Comparison = "Greater_or_Eq"} else {Ratio_Comparison = "Less"}
        #add to list---could recode this to avoid having this extra step and list...
        list_of_ratio_comparisons = append(list_of_ratio_comparisons,Ratio_Comparison)
      
        }
      }
      #calculate pvals
      #Divide # of greater/eq comparisons by P-total number of perm fails--->instances where we had no inter coil pairs
      greater_eq_div_total = (length(grep("Greater_or_Eq",list_of_ratio_comparisons)))/(P-length(list_of_perm_fails))
      
      #add to dataframe
      message(paste0("Done: ", Structure_1))
      all = cbind(Structure_1, greater_eq_div_total)
      permtest_inter_intra_coil_type_Distortion_Corrected_pvals_means= rbind(permtest_inter_intra_coil_type_Distortion_Corrected_pvals_means,all)
    }
      
