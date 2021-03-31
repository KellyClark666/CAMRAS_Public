library(neurobase)
library(ANTsR)
library(extrantsr)
library(fslr)
library(parallel)
library(WhiteStripe)



id.list = list.files(pattern = "-",recursive = FALSE, full.names = FALSE)

do.reg = function(id){
  
  sites_all = list.files(path = id, recursive = FALSE, full.names = TRUE)
  sites_no_hop = sites_all[which(!grepl("Hopkins",sites_all))] 
  sites_hop = sites_all[which(grepl("Hopkins",sites_all))] 
  
  for (i in 1:length(sites_no_hop)){
 
  t1 = system(paste("find ", paste0(sites_no_hop[i],"/analysis")," -name MPRAGE_SAG_TFL*n4.nii.gz -o -name *3D_T1_MPRAGE_CAVS*n4.nii.gz"), intern = TRUE)
  t1n4_scan1=t1[which(!grepl("*a_n4.nii.gz",t1))]
  t1n4_scan2=t1[which(grepl("*a_n4.nii.gz",t1))]
  t1n4_scan1_ND=t1n4_scan1[which(grepl("*ND",t1n4_scan1))]
  t1n4_scan1_corr = t1n4_scan1[which(!grepl("*ND",t1n4_scan1))]
  t1n4_scan2_ND=t1n4_scan2[which(grepl("*ND",t1n4_scan2))]
  t1n4_scan2_corr = t1n4_scan2[which(!grepl("*ND",t1n4_scan2))]
  
  flair = system(paste0("find ",paste0(sites_no_hop[i],"/analysis")," -name *FLAIR_SAG_VFL*n4.nii.gz -o -name *FLAIR_SAG-VFL*n4.nii.gz "), intern = TRUE)
  flairsn4_scan1_all=flair[which(!grepl("*a_n4.nii.gz",flair))]
  flairsn4_scan2_all=flair[which(grepl("*a_n4.nii.gz",flair))]
  flairsn4_scan1_ND=flairsn4_scan1_all[which(grepl("*ND",flairsn4_scan1_all))]
  flairsn4_scan1_corr = flairsn4_scan1_all[which(!grepl("*ND",flairsn4_scan1_all))]
  flairsn4_scan2_ND=flairsn4_scan2_all[which(grepl("*ND",flairsn4_scan2_all))]
  flairsn4_scan2_corr = flairsn4_scan2_all[which(!grepl("*ND",flairsn4_scan2_all))]

  flair1_to_mprage1 = antsRegistration(fixed = oro2ants(readnii(t1n4_scan1_corr)), 
                                       moving = oro2ants(readnii(flairsn4_scan1_corr)),
                                       typeofTransform = "Rigid")
  flair1_reg = antsApplyTransforms(fixed=oro2ants(readnii(t1n4_scan1_corr)), 
                                   moving=oro2ants(readnii(flairsn4_scan1_corr)),
                                  transformlist = flair1_to_mprage1$fwdtransforms,
                                  interpolator ="WelchWindowedSinc")
  #antsImageWrite(flair1_reg, paste0(dirname(t1n4_scan1_corr),"/flair1_reg_to_mprage1.nii.gz"))
  
  flair2_to_mprage2 = antsRegistration(fixed = oro2ants(readnii(t1n4_scan2_corr)), 
                                       moving = oro2ants(readnii(flairsn4_scan2_corr)),
                                       typeofTransform = "Rigid")
  flair2_reg = antsApplyTransforms(fixed=oro2ants(readnii(t1n4_scan2_corr)), moving=oro2ants(readnii(flairsn4_scan2_corr)),
                                   transformlist = flair2_to_mprage2$fwdtransforms,
                                   interpolator ="WelchWindowedSinc")
  #antsImageWrite(flair2_reg, paste0(dirname(t1n4_scan2_corr),"/flair2_reg_to_mprage2.nii.gz"))
  
  flair1ND_to_mprage1ND = antsRegistration(fixed = oro2ants(readnii(t1n4_scan1_ND)), 
                                       moving = oro2ants(readnii(flairsn4_scan1_ND)),
                                       typeofTransform = "Rigid")
  flair1ND_reg = antsApplyTransforms(fixed=oro2ants(readnii(t1n4_scan1_ND)), 
                                   moving=oro2ants(readnii(flairsn4_scan1_ND)),
                                   transformlist = flair1ND_to_mprage1ND$fwdtransforms,
                                   interpolator ="WelchWindowedSinc")
  #antsImageWrite(flair1ND_reg, paste0(dirname(t1n4_scan1_ND),"/flair1ND_reg_to_mprage1ND.nii.gz"))
  
  flair2ND_to_mprage2ND = antsRegistration(fixed = oro2ants(readnii(t1n4_scan2_ND)), 
                                           moving = oro2ants(readnii(flairsn4_scan2_ND)),
                                           typeofTransform = "Rigid")
  flair2ND_reg = antsApplyTransforms(fixed=oro2ants(readnii(t1n4_scan2_ND)), 
                                     moving=oro2ants(readnii(flairsn4_scan2_ND)),
                                     transformlist = flair2ND_to_mprage2ND$fwdtransforms,
                                     interpolator ="WelchWindowedSinc")
  #antsImageWrite(flair2ND_reg, paste0(dirname(t1n4_scan2_ND),"/flair2ND_reg_to_mprage2ND.nii.gz"))
  }
  
  for (i in 1:length(sites_hop)){
    
    t1 = system(paste("find ", paste0(sites_no_hop[i],"/analysis")," -name MPRAGE_SAG_TFL*n4.nii.gz -o -name *3D_T1_MPRAGE_CAVS*n4.nii.gz"), intern = TRUE)
    t1n4_scan1=t1[which(!grepl("*a_n4.nii.gz",t1))]
    t1n4_scan2=t1[which(grepl("*a_n4.nii.gz",t1))]
    t1n4_scan1_corr = t1n4_scan1[which(!grepl("*ND",t1n4_scan1))]
    t1n4_scan2_corr = t1n4_scan2[which(!grepl("*ND",t1n4_scan2))]
   
    flair = system(paste0("find ",paste0(sites_no_hop[i],"/analysis")," -name *FLAIR_SAG_VFL*n4.nii.gz -o -name *FLAIR_SAG-VFL*n4.nii.gz "), intern = TRUE)
    flairsn4_scan1_all=flair[which(!grepl("*a_n4.nii.gz",flair))]
    flairsn4_scan2_all=flair[which(grepl("*a_n4.nii.gz",flair))]
    flairsn4_scan1_corr = flairsn4_scan1_all[which(!grepl("*ND",flairsn4_scan1_all))]
    flairsn4_scan2_corr = flairsn4_scan2_all[which(!grepl("*ND",flairsn4_scan2_all))]
    
    flair1_to_mprage1 = antsRegistration(fixed = oro2ants(readnii(t1n4_scan1_corr)), 
                                         moving = oro2ants(readnii(flairsn4_scan1_corr)),
                                         typeofTransform = "Rigid")
    flair1_reg = antsApplyTransforms(fixed=oro2ants(readnii(t1n4_scan1_corr)), 
                                     moving=oro2ants(readnii(flairsn4_scan1_corr)),
                                     transformlist = flair1_to_mprage1$fwdtransforms,
                                     interpolator ="WelchWindowedSinc")
    antsImageWrite(flair1_reg, paste0(dirname(t1n4_scan1_corr),"/flair1_reg_to_mprage1.nii.gz"))
    
    flair2_to_mprage2 = antsRegistration(fixed = oro2ants(readnii(t1n4_scan2_corr)), 
                                         moving = oro2ants(readnii(flairsn4_scan2_corr)),
                                         typeofTransform = "Rigid")
    flair2_reg = antsApplyTransforms(fixed=oro2ants(readnii(t1n4_scan2_corr)), moving=oro2ants(readnii(flairsn4_scan2_corr)),
                                     transformlist = flair2_to_mprage2$fwdtransforms,
                                     interpolator ="WelchWindowedSinc")
    antsImageWrite(flair2_reg, paste0(dirname(t1n4_scan2_corr),"/flair2_reg_to_mprage2.nii.gz"))
  }  
  
}

mclapply(id.list, do.reg, mc.cores = length(id.list))
