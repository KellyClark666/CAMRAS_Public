library(fslr)
library(neurobase)
library(ANTsR)
library(extrantsr)
library(scales)
library(parallel)


dirs= system(paste("find ./Data/ -name mass"),intern=TRUE) 

do.first=function(j){

  #list all mass skull stripped images
  imgs = list.files(path = dirs[j], pattern = "*_brain.nii.gz",full.names = TRUE)
  
  #subest images to create list of only the mass skull stripped mprages
  t1s = imgs[which(grepl("MPRAGE",imgs))]
  
  for (i in 1:length(t1s)){
    #create FIRST output dir
    first.dir = paste0(dirname(dirs[j]),"/FIRST")
    
    #Create output dir in FIRST dir for each of the 4 images (scan/rescan and ND scan/rescan)
    output.dir = paste0(basename(gsub("_brain.nii.gz","",(t1s[i]))))
    output.dir_path = file.path(first.dir, output.dir )
    dir.create(output.dir_path)
    
    print(paste0("Directory created for: ", t1s[i]))
    
    #R wrapper for FIRST command 
    system(paste("run_first_all -b -s R_Thal,L_Thal -i",t1s[i],"-o", paste0(output.dir_path,"/",basename(nii.stub(t1s[i])),"_thalamus")))
  }
}
j = 1:length(dirs)
mclapply(j, do.first, mc.cores = 50)
