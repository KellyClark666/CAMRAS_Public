###Mean(or Median,adjust accordingly) inter/intra coil ratios from Observed Data###    
#this is very ugly, need to tidy this up and get rid of some for loops 

library(tidyverse)
library(reshape2)
library(dplyr)

 #BWH
 subjects_1 = c("Subject","list","for","BWH")
 #Penn
 subjects_2 = c("Subject","list","for","Penn")

df = read.csv("./path/to/df/containing/volumetric info")
df = melt(df)

inter_intra_coil_type_actual_ratios = data.frame()
for(k in 1:length(struc_list)){
  Structure = struc_list[k]  
  struc_df = BWH_Penn_D %>% filter(variable==Structure)
  coil_intra_list = vector()
  coil_inter_list = vector()
  for (j in 1:length(site_list)){
    Site1 = site_list[j]
    vols_df = struc_df %>% filter(Site==Site1)
    
  if (Site1=="Penn"){
    print("ok")
    for (i in 1:length(subjects_2)){
      Subject = subjects_2[i]
      
      pair = vols_df %>% filter(Subj==Subject)
      if (pair[,4,1][1]==pair[,4,1][2]){
        intracoil_diffs = dist(pair$value)
        message(paste0("Intra-", Subject,":",intracoil_diffs))
        coil_intra_list = append(coil_intra_list, intracoil_diffs)
        
      }else{
        intercoil_diffs = dist(pair$value)
        message(paste0("Inter-",Subject,":",intercoil_diffs))
        
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
    avg_coil_intra_both_sites = mean(coil_intra_list,na.rm=TRUE)
    avg_coil_inter_both_sites = mean(coil_inter_list,na.rm=TRUE)
    ratio_inter_intra_avg_coil_diffs = avg_coil_inter_both_sites/avg_coil_intra_both_sites
    all = cbind(Structure,avg_coil_intra_both_sites,avg_coil_inter_both_sites,ratio_inter_intra_avg_coil_diffs)
    inter_intra_coil_type_actual_ratios = rbind(inter_intra_coil_type_actual_ratios,all)    
  }
  
  
  ###Mean inter/intra site ratios from Observed Data (adjust D and ND / including and excluding hopkins accordingly###  


#read in data, melt down to long form
cam = read.csv("/path/to/data")

#filter/subset to get Distortion and non-dist df's 
cam_1  = cam %>% filter(Type=="D" | Site=="Hopkins")
cam_2 = cam %>% filter(Type=="ND")
data_D_Hop = melt(cam_1)
data_ND_Hop = melt(cam_2)

#Subjects
subjects_1 = c("list","of","subjects")

#list of structures
struc_list = c("JLF.Thal.Vol","First.Thal.Vol","Mimosa.Lesion.Vol",
               "JLF.White.Matter.Vol","JLF.Grey.Matter.Vol","FAST.TBV","Atropos.Grey.Matter.Vol")

#include/exclude hopkins

D_and_Hopkins_df = data_D_Hop #this is very pointless...IDK why I did this
D_df = data_D_Hop %>% filter(Type=="D")

ND_incl_hopkins = data_ND_Hop
ND_excl_hopkins = data_ND_Hop %>% filter(Site!="Hopkins")

set.seed(1979)

#Initialize empty DF to store results
ND_and_Hop_inter_intra_ratio=data.frame()

#loop thru structures (we'll end up with a p-value for each structure)
for (k in 1:length(struc_list)){  
  Structure = struc_list[k] 
  
  #Subset original df based on "Structure"
  vols_df = ND_incl_hopkins %>% filter(variable==Structure)
  #get rid of any NAs
  vols_df = vols_df[order(vols_df$Site),]
  across_subj_avg_inter = vector()
  across_subj_avg_intra = vector()
  
    #loop through each Subject within each permutation (so I can look outputs for Perm[,l] for each subject)
    for(i in 1:length(subject_list)){
      Subject = subject_list[i]
      #if subject =="03-002" | "03-001" do things differently 
      if (Subject=="03-002"){
        
        #subset vols_df by "Subject"
        new_df = vols_df %>% filter(Subj==Subject)
        
        #create vector of all possible inter/intra site abs differences (done using dist functin)
        diff_vector = as.numeric(dist(new_df$value))  
        
        #subset vector into intersite and intrasite pairs only (based on vols_df site labels ordered alphabetically)
        abs_diff_for_each_intersite_pair = diff_vector[c(2:9,11:14)]
        abs_diff_for_each_intrasite_pair = diff_vector[-c(2:9,11:14)]
        
      } else if (Subject=="03-001"){
        new_df = vols_df %>% filter(Subj==Subject)

        diff_vector = as.numeric(dist(new_df$value))  
        
        abs_diff_for_each_intersite_pair = diff_vector[c(2:5)] 
        abs_diff_for_each_intrasite_pair = diff_vector[-c(2:5)] 
        
      } else {
        new_df = vols_df %>% filter(Subj==Subject)

        diff_vector = as.numeric(dist(new_df$value)) 
        
        abs_diff_for_each_intersite_pair = diff_vector[c(2:13,15:22,24:27)]
        abs_diff_for_each_intrasite_pair = diff_vector[-c(2:13,15:22,24:27)]
        
      }
      
      #average across intersite and intrasite difference within subject 
      average_abs_diff_across_intersite_pairs = mean(abs_diff_for_each_intersite_pair)
      average_abs_diff_across_intrasite_pairs = mean(abs_diff_for_each_intrasite_pair)
      
      #add within subject avg inter/intra diffs to a list
      #you'll end up with a list of avg inter/intra diffs per subjects for whatever permutation you're on (Permutation[,l])
      across_subj_avg_inter = append(across_subj_avg_inter,average_abs_diff_across_intersite_pairs)
      across_subj_avg_intra = append(across_subj_avg_intra,average_abs_diff_across_intrasite_pairs)
    }
    
    #After you get that list ^^ calculate inter/intra ratio by averaging those avg inter/intra ^ across all subjects then dividing that inter/intra avg
   Mean_ABS_Diff_InterSite = mean(across_subj_avg_inter)
   Mean_ABS_Diff_IntraSite = mean(across_subj_avg_intra)
  ratio = Mean_ABS_Diff_InterSite/Mean_ABS_Diff_IntraSite
  
  #add to dataframe
  all = cbind(Structure, Mean_ABS_Diff_InterSite, Mean_ABS_Diff_IntraSite, ratio)
  ND_and_Hop_inter_intra_ratio= rbind(ND_and_Hop_inter_intra_ratio,all)
  
}
