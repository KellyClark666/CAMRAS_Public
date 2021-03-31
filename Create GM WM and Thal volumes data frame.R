library(neurobase)
library(vroom)
library(parallel)


dirs= system(paste("find . -name analysis"),intern=TRUE)  

do.vol = function(i){
  
  gmwm_masks = list.files(path = dirs[i],pattern = "jlf_wmgm_mask.nii.gz", 
                     recursive=TRUE, full.names=TRUE)
  thalamus_masks = list.files(path = dirs[i],pattern = "thalamus_all_none_firstseg.nii.gz",
                              recursive = TRUE, full.names = TRUE)
  
 
  df = data.frame()
  for (j in 1:length(gmwm_masks)){
    labels =  strsplit(dirname(dirs[i]),"/")[[1]]
    Subj = labels[2]
    Site = labels[3]
    labels2 = strsplit(gmwm_masks[j],"/")[[1]]
    Scan =  labels2[6]
    gmwm_mask = readnii(gmwm_masks[j])
    grey_antslab = table(gmwm_mask[gmwm_mask == 2])
    white_antslab = table(gmwm_mask[gmwm_mask==1])
    vreslab = voxres(gmwm_mask, units = "cm")
    GreyMatter_Vol = grey_antslab * vreslab
    WhiteMatter_Vol = white_antslab * vreslab
    
    thal = readnii(thalamus_masks[j])
    thal_lab_1 = (table(thal[thal==10]))
    thal_lab_2 = (table(thal[thal==49]))
    thal_lab = voxres(thal, units ="cm")
    Thal_Vol = (thal_lab*thal_lab_1)+(thal_lab*thal_lab_2)
            
    
    all_together = cbind(Subj,Site,Scan,GreyMatter_Vol, WhiteMatter_Vol, Thal_Vol)
    df = rbind(df, all_together)
    rownames(df)=c()
    write.csv(df,paste0(dirs[i],"/volumes_check.csv"))
  }
  
  }
i = 1:41
mclapply(i,do.vol,mc.cores = 40)

Volumes_Files=list.files(pattern="volumes.csv",recursive=TRUE,full.names=TRUE)

All_Volumes = vroom(Volumes_Files)
write.csv(All_Volumes, "GMWM_Thal_Volumes.csv")
