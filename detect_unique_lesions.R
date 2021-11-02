#Kelly Clark
#11/01/2021

#This code goes through each connected component (lesion) in the sum of all segmentation images for a given subject
#and checks to see if the lesion exists in each individual segmentation
#if the lesion is disappearing (ie., does not exist or was not segmented)
#the visit where the disappearing lesion occurred and the index of a voxel where this exists in the sum image
#is recorded in a df to be checked visually for errors or actual cases of "disappearing lesions"

library(neurobase)
library(ANTsR)
library(extrantsr)
library(fslr)
library(parallel)
setwd("/path/to/manual/segmentations")


# list all of the manual segmentation images that have CC's labelled
scans = list.files("./registered_FLAIRS_and_FLAIR_segmentations",pattern="_lesion_binary_mask_registered.nii.gz",full.names = TRUE)
subject = c("01-001","01-002","01-003","02-001","02-002","02-003","03-001","03-002","04-001","04-002","04-003")

do.lesion = function(i){
  #Extract ID and Subj for print statements
  label = gsub("_FLAIR_SAG_VFL_lesion_binary_mask_registered.nii.gz","",strsplit(scans[i],"/")[[1]][3])
  subj = strsplit(label,"_")[[1]][1]
  Site_Visit = strsplit(label,"_")[[1]][2]
  #list all of the "sum" files (sum of all manual segmentations for given subject)
  sum_path = paste0("./sum_registered_lesion_masks/cc_sum_lesion_masks/",subj,"_lesion_masks_sum_CC.nii.gz")
  sum = readnii(sum_path)
  img = readnii(scans[i])

  #get total number of connnected components
  number_of_CC = 1:max(sum)

  #loop through each connected component (lesion)
  test = list() 
  
   for (l  in number_of_CC){

  #subset individual lesion of interest by CC label  
  lesion = which(sum==l, arr.ind=TRUE)

  #loop through each voxel index within lesion 
  for (j in 1:nrow(lesion)){

  index = lesion[j,]
 
  x = (index[1])
  y = (index[2])
  z = (index[3])
  
  #print messages to display if condition is met
  #unique = paste0("Lesion ", l, " is unique")
  not_unique = paste0("Lesion ", l, " from ",subj, " is not unique")
   
  #if any of the voxels in the lesion of interest correspond to an 
  #intensity of 8 (or 4 or 6 for 03-001/03-002) in the sum of all the segmentation masks
  #then lesion is not unique, stop and current iteration and move to next lesion
  if (img[x,y,z]==1){
      print(not_unique)
    test = c(test,l)
    not_unique_list=as.numeric(unlist(test))
      stop=TRUE
       break}
  }
  
  unique_list = number_of_CC[!(number_of_CC %in% test)]
  unique_lesion_df = data.frame()
  for (z in 1:length(unique_list)){
    lesion_label = unique_list[z]
    lesion = which(sum==unique_list[z], arr.ind=TRUE)
    index = lesion[1,]
    x = (index[1])
    y = (index[2])
    z = (index[3])
    ok = cbind(subj,lesion_label,x,y,z)
    unique_lesion_df = rbind(unique_lesion_df,ok)
     unique_lesion_df[!duplicated(unique_lesion_df[c(1,2)]),]
  }
            
  
   }
  rownames(unique_lesion_df)=NULL
  # write.csv(unique_lesion_df,paste0("./",label,"_unique_lesion_df.csv"))
            
          
}

i = 1:length(scans)
mclapply(i, do.lesion,mc.cores=11)


