####Permutation Test: Comparison of ratios of mean inter/mean intra site diffs after permutation to observed data ratio 
####Kelly Clark 
####July 2, 2021
####just a bunch of for loops LOL 

#read in data, melt down to long form
cam_1 = read.csv("/path/to/data")
data = melt(cam_1)

#Subjects
subjects_1 = c("list","of","subject","ids")

struc_list = c("list","of","structures","we","have","data","for")

set.seed(1979)

#Initialize empty DF to store results
perm_inter_intra_pvals=data.frame()

#loop thru structures (we'll end up with a p-value for each structure)
for (k in 1:length(struc_list)){  
  Structure = struc_list[k] 
  
  #Subset original df based on "Structure"
  vols_df = data %>% filter(variable==Structure)
  #get rid of any NAs
  vols_df = vols_df[order(vols_df$Site),]
  
  #get inter/intra ratio for Structure from actual data
  inter_intra_ratio_actual_data_df = actual_data_avg_of_abs_diffs_inter_intra_ratio %>% filter(Structure1==Structure)
  inter_intra_actual_ratio = inter_intra_ratio_actual_data_df$Incl.Hopkins.Ratio
  
  #number of permutations 
  P = 10000
  
  #initialize a bunch of empty lists to store ouputs from each iteration of for loop (lol!)
  list_of_ratios = vector()
  list_of_ratio_comparisons = vector()
  
  for(l in 1:P){
    #more empty lists :')
    across_subj_avg_inter = vector()
    across_subj_avg_intra = vector()
    
    #loop through each Subject within each permutation (so I can look outputs for Perm[,l] for each subject)
    for(i in 1:length(subject_list)){
      Subject = subject_list[i]
      #if subject =="03-002" | "03-001" do things differently 
      if (Subject=="03-002"){
        
        #subset vols_df by "Subject"
        new_df = vols_df %>% filter(Subj==Subject)
        
        #specify sample size (however many imgs/measurements were acquired for Subject)
        n = length(new_df$Site)
        
        #specify variable to shuffle (volumes)
        variable = new_df$value
        
        #initialize empty matrix...I'm not sure this is necessary since were looping through subject within the permutation for loop
        PermSamples = matrix(0, nrow=n, ncol=P)
        PermSamples[,l] = sample(variable, size=n, replace=FALSE)
        
        #create vector of differences for all possible pairs of images
        PermDF = as.data.frame(PermSamples[,l])
        
        #create vector of all possible inter/intra site abs differences (done using dist functin)
        diff_vector = as.numeric(dist(PermDF$`PermSamples[, l]`))  
        
        #subset vector into intersite and intrasite pairs only (based on vols_df site labels ordered alphabetically)
        abs_diff_for_each_intersite_pair = diff_vector[c(2:9,11:14)]
        abs_diff_for_each_intrasite_pair = diff_vector[-c(2:9,11:14)]
        
      } else if (Subject=="03-001"){
        new_df = vols_df %>% filter(Subj==Subject)
        n = length(new_df$Site)
        variable = new_df$value
        
        PermSamples = matrix(0, nrow=n, ncol=P)
        PermSamples[,l] = sample(variable, size=n, replace=FALSE)
        
        PermDF = as.data.frame(PermSamples[,l])
        diff_vector = as.numeric(dist(PermDF$`PermSamples[, l]`))  
        
        abs_diff_for_each_intersite_pair = diff_vector[c(2:5)] 
        abs_diff_for_each_intrasite_pair = diff_vector[-c(2:5)] 
        
      } else {
        new_df = vols_df %>% filter(Subj==Subject)
        n = length(new_df$Site)
        variable = new_df$value
        
        PermSamples = matrix(0, nrow=n, ncol=P)
        PermSamples[,l] = sample(variable, size=n, replace=FALSE)
        
        PermDF = as.data.frame(PermSamples[,l])
        diff_vector = as.numeric(dist(PermDF$`PermSamples[, l]`)) 
        
        abs_diff_for_each_intersite_pair = diff_vector[c(2:13,15:22,24:27)]
        abs_diff_for_each_intrasite_pair = diff_vector[-c(2:13,15:22,24:27)]
        
      }
      
      #average across intersite and intrasite difference within subject, I don  
      average_abs_diff_across_intersite_pairs = mean(abs_diff_for_each_intersite_pair)
      average_abs_diff_across_intrasite_pairs = mean(abs_diff_for_each_intrasite_pair)
      
      #add within subject avg inter/intra diffs to a list
      #you'll end up with a list of avg inter/intra diffs per subjects for whatever permutation you're on (Permutation[,l])
      across_subj_avg_inter = append(across_subj_avg_inter,average_abs_diff_across_intersite_pairs)
      across_subj_avg_intra = append(across_subj_avg_intra,average_abs_diff_across_intrasite_pairs)
    }
    
    #After you get that list ^^ calculate inter/intra ratio by averaging those avg inter/intra ^ across all subjects then dividing that inter/intra avg
    ratio = mean(across_subj_avg_inter)/mean(across_subj_avg_intra)
    #add to a list of ratios
    list_of_ratios = append(list_of_ratios,ratio)
    #is ratio after permutation >= observed ratio? or not?
    if (ratio>=inter_intra_actual_ratio){ 
      Ratio_Comparison = "Greater_or_Eq"} else {Ratio_Comparison = "Less"}
    
    #add to list
    list_of_ratio_comparisons = append(list_of_ratio_comparisons,Ratio_Comparison)
  }
  #calculate p value for "structure[k]"Structure" by dividing # of ratios that are >= observed ratio by total ratios (10000 bc P =10000)
  greater_eq_div_total = (length(grep("Greater_or_Eq",list_of_ratio_comparisons)))/P
  
  #add to dataframe
  all = cbind(Structure, greater_eq_div_total)
  perm_inter_intra_pvals= rbind(perm_inter_intra_pvals,all)
  
}







