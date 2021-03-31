#Kelly Clark
#March 30, 2020
#N4 processing of CAMRAS data

library(neurobase)
library(ANTsR)
library(extrantsr)
library(fslr)
library(parallel)
library(WhiteStripe)



id.list = list.files(pattern = "-",recursive = FALSE, full.names = FALSE)
do.n4 = function(id){
t1s = system("find  -name MPRAGE_SAG_TFL*.nii.gz -o -name *3D_T1_MPRAGE_CAVS*.nii.gz", intern = TRUE)
t1s = t1s[which(grepl("/NIFTI*",t1s))]

flairs = system("find . -name *FLAIR_SAG_VFL*.nii.gz -o -name *FLAIR_SAG-VFL*.nii.gz ", intern = TRUE)
flairs = flairs[which(grepl("/NIFTI*",flairs))]

all.imgs = c(t1s, flairs)

for (i in 1:length(all.imgs)){
    n4 = bias_correct(all.imgs[i], correction = "N4")
    subj_dir = dirname(dirname(all.imgs[i]))
    dir.create(file.path(subj_dir,"analyis"))
    writenii(n4, paste0(subj_dir,"/analysis/", nii.stub(basename( all.imgs[i])),"_n4.nii.gz"))
    print("Done")
  }
}

mclapply(id.list, do.n4, mc.cores = length(id.list))
