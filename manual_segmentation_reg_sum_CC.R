##Kelly Clark
##10/28/21

##This code creates a binary mask from lesion maps,
##registers all binary masks to a common image within each subject
##Sums the segmentation masks into one image for each subject
##And labels connected components for each summed segmentation image 

library(neurobase)
library(ANTsR)
library(extrantsr)
library(fslr)
library(parallel)
setwd("/path/to/manual/segmentations")

##create binary mask from lesion maps
seg_maps = list.files(".",pattern = "_lesion_mask.nii.gz",full.names = TRUE,recursive = TRUE )
for (i in 2:length(seg_maps)){
  
  seg_map = readnii(seg_maps[i])
  mask = seg_map>0
  writenii(mask,paste0(gsub("_mask.nii.gz","_binary_mask.nii.gz",seg_maps[i])))
  message(paste("Done: ",seg_maps[i]))
  
}

####Register all segmentations to each subjects's BWH scan 1

#Subject list to loop through
subject = c("01-001","01-002","01-003","02-001","02-002","02-003","03-001","03-002","04-001","04-002","04-003")

#save all segmentations to a list
mprage_manual_seg_masks = list.files(".",pattern = "MPRAGE_SAG_TFL_lesion_binary_mask.nii.gz",recursive=TRUE,full.names=TRUE)
mprage_manual_seg_1 = mprage_manual_seg_masks[which(!grepl("a/",mprage_manual_seg_masks))]
mprage_manual_seg_2 = mprage_manual_seg_masks[which(grepl("a/",mprage_manual_seg_masks))]

flair_manual_seg_masks = list.files(".",pattern = "FLAIR_SAG_VFL_lesion_binary_mask.nii.gz",recursive=TRUE,full.names=TRUE)
flair_manual_seg_1 = flair_manual_seg_masks[which(!grepl("a/",flair_manual_seg_masks))]
flair_manual_seg_2 = flair_manual_seg_masks[which(grepl("a/",flair_manual_seg_masks))]

do.stuff = function(j){
  message(paste0("Looking for masks: ",subject[j])) 
  
  manual_1_all = flair_manual_seg_1[which(grepl(subject[j],flair_manual_seg_1))]
  manual_2_all = flair_manual_seg_2[which(grepl(subject[j],flair_manual_seg_2))]
  manual_1 = manual_1_all[which(!grepl(paste0(subject[j],"_BWH"),manual_1_all))]
  manual_fixed = manual_1_all[which(grepl(paste0(subject[j],"_BWH"),manual_1_all))]
  
  #list of manual segmentations to be registeres
  manual = c(manual_2_all,manual_1)
  
  #BWH Scan 1 segmentation doesn't need to be registered
  fixed = manual_1_all[which(grepl("_BWH",manual_1_all))]
  message(paste0("Masks and images found: ",subject[j])) 
  
  #Scans to be registered to BWH Scan 1 FLAIR
  moving_imgs_paths = paste0(dirname(manual),"/FLAIR_SAG_VFL.nii.gz")
  
  #BWH Scan 1 FLAIR image path (what we are registering everything to)
  fixed_img_path = paste0(dirname(fixed),"/FLAIR_SAG_VFL.nii.gz")
  
  message(paste0("For Loop Time ye ye: ",subject[j])) 
  
  for (i in 1:length(moving_imgs_paths)){

    #register all flair images for subject to BWH scan 1 flair
    reg_to_BWH = antsRegistration(fixed = oro2ants(readnii(fixed_img_path)),
                                  moving = oro2ants(readnii(moving_imgs_paths[i])),
                                  typeofTransform = "Rigid")
    #apply transform above to flair images to actually make the image be in the same space as BWH scan1
    reg_flair = antsApplyTransforms(fixed=oro2ants(readnii(fixed_img_path)),
                                    moving=oro2ants(readnii(moving_imgs_paths[i])),
                                    transformlist = reg_to_BWH$fwdtransforms,
                                    interpolator ="WelchWindowedSinc")
    message(paste0("Registered: ",moving_imgs_paths[i]))

    #register manual segmentations to BWH scan1 by applying transform from image registration
    reg_manual_segmentation = antsApplyTransforms(fixed=oro2ants(readnii(fixed_img_path)),
                                                  moving=oro2ants(readnii(manual[i])),
                                                  transformlist = reg_to_BWH$fwdtransforms,
                                                  interpolator ="nearestNeighbor")

    #pull subject, site and session label from file path to write out images with
    labels = strsplit(moving_imgs_paths[i],"/")[[1]]
    x = labels[2]

    #write out registered flair images and registered segmentations
    antsImageWrite(reg_flair, paste0("./registered_FLAIRS_and_FLAIR_segmentations/",x,"_FLAIR_SAG_VFL_registered.nii.gz"))
    antsImageWrite(reg_manual_segmentation, paste0("./registered_FLAIRS_and_FLAIR_segmentations/",x,"_FLAIR_SAG_VFL_lesion_binary_mask_registered.nii.gz"))

  }
  writenii(readnii(fixed_img_path),paste0("./registered_FLAIRS_and_FLAIR_segmentations/",subject[j],"_BWH_FLAIR_SAG_VFL.nii.gz"))
  writenii(readnii(manual_fixed),paste0("./registered_FLAIRS_and_FLAIR_segmentations/",subject[j],"_BWH_FLAIR_SAG_VFL_lesion_binary_mask_registered.nii.gz"))
  message(paste0("Done: ",subject[j])) 
}
j= 2:length(subject)
mclapply(j,do.stuff,mc.cores=11)

####Sum segmentations for each subject
flair_manual_seg_masks = list.files(".",pattern = "FLAIR_SAG_VFL_lesion_binary_mask_registered.nii.gz",recursive=TRUE,full.names=TRUE)
do.sum = function(j){
  manual_all = flair_manual_seg_masks[which(grepl(subject[j],flair_manual_seg_masks))]
  
  img = 0 
    for (i in 1:length(manual_all)){
     img = readnii(manual_all[i]) + img
    }
  writenii(img,paste0("./sum_registered_lesion_masks/",subject[j],"_lesion_masks_sum.nii.gz"))
  
}
j = 2:length(subject)
mclapply(j,do.sum,mc.cores=11)

####CC label each image
#list out all registered manual segmentation leson masks
flair_manual_seg_masks = list.files("./sum_registered_lesion_masks/",pattern = "_lesion_masks_sum.nii.gz",recursive=TRUE,full.names=TRUE)
do.cc = function(j){
  
  #read in lesion mask
  flair_manual_seg_mask = readnii(flair_manual_seg_masks[j])
  message(paste("Lesion Mask Read In: ",flair_manual_seg_masks[j]))
  #label connected components 
  cc_labelled_lesion_mask = label_mask(flair_manual_seg_mask)
  message(paste("CC labelling finished: ",flair_manual_seg_masks[j]))
  #write out cc labelled image
  writenii(cc_labelled_lesion_mask,paste0("./sum_registered_lesion_masks/cc_sum_lesion_masks/",
                                          nii.stub(basename(flair_manual_seg_masks[j])),"_CC.nii.gz"))
  message(paste("CC Labelled Image Written Out: ",flair_manual_seg_masks[j]))
  
  }
j = 1:length(flair_manual_seg_masks)
mclapply(j,do.cc,mc.cores=50)

