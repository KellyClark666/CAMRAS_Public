library(neurobase)
library(parallel)
library(vroom)




dirs = system(paste("find . -name FAST"),intern=TRUE) 

do.tbv = function(i){

#list all output files from FAST
FAST = list.files(path = paste0(dirs[i]), recursive=TRUE,full.names=TRUE)

#subset list of output files 
FAST_TB_seg = FAST[which(grepl("*_brain_seg.nii.gz",FAST))]

df = data.frame()

for (j in 1:length(FAST_TB_seg)){
  labels=strsplit(FAST_TB_seg[j],"/")[[1]]
  Subj = labels[2]
  Site = labels[3]
  Scan = gsub("_brain_seg.nii.gz","",labels[6])  
  seg_img = readnii(FAST_TB_seg[j])
  
  #Creat binary mask. CSF intensity = 1, all other tissue is > 1
  Total_Brain = seg_img>1
  Total_Brain_Lab = table(Total_Brain[Total_Brain==1])
  vreslab = voxres(seg_img,units = "cm")
  TBV = Total_Brain_Lab*vreslab
  
  all_together = cbind(Subj,Site,Scan,TBV)
  df = rbind(df, all_together)
  #rownames(df)=c()
  write.csv(df,paste0(dirname(dirs[i]),"/TBV.csv"))
  
  
}
#print(df)
}

i = 1:length(dirs)
do.tbv(i)

files = list.files(pattern = "TBV.csv", recursive = TRUE, full.names = TRUE)
combined_df = vroom(files)
write.csv(combined_df,"./TBV_all_subjects.csv")
