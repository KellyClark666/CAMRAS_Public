library(ANTsR)
library(extrantsr)
library(neurobase)
library(parallel)


dirs= system(paste("find ./Data/ -name mass"),intern=TRUE) 

do.jlf = function(x) {
  #list all mass skull stripped images
  imgs = list.files(path = dirs[x], pattern = "*_brain.nii.gz",
                    recursive = TRUE, full.names = TRUE)
  
  #subset imgs to get mass skull stripped mprages only
  t1s = imgs[which(grepl("MPRAGE",imgs))]
  
  for (i in 1:length(t1s)){
 
    #Create JLF ouput dirs
    JLFdir = paste0(dirname(dirname(t1s[i])),"/JLF_WMGM")
    outdir_regatlas = paste0(gsub("_brain.nii.gz","",(basename(t1s[i]))),"/muse_to_t1")
    outdir_regwmgm = paste0(gsub("_brain.nii.gz","",(basename(t1s[i]))),"/muse_wmgm_to_t1")
    output.dir_path_atlas = file.path(JLFdir, outdir_regatlas )
    output.dir_path_wmgm = file.path(JLFdir, outdir_regwmgm)
    dir.create(output.dir_path_atlas, recursive = TRUE)
    dir.create(output.dir_path_wmgm, recursive = TRUE)
    print(paste0("Dirs created for ", t1s[i]))
    

    for(j in 1:10){
      
      t1 = readnii(t1s[i])
      
      #full brain atlases
      muse = paste0("/project/Melissa_stuff/MUSE_Templates/WithCere/Template", j, ".nii.gz")
      
      #WM/GM atlases 
      wmgm = paste0("/project/Melissa_stuff/MUSE_Templates/WithCere/Template", j, "_label_WMGM.nii.gz")
      
      # Register muse full brain atlases to camras images
      print(paste0("Registering muse atlas to ", t1s[i]))
      atlas_to_image = registration(filename = muse,
                                    template.file = t1,
                                    typeofTransform = "SyN", remove.warp = FALSE)
      
      # Apply registration to wm/gm atlases and full brain atlases 
      print(paste0("Applying transforms to brains and wmgm for ", t1s[i]))
      muse_reg = antsApplyTransforms(fixed = oro2ants(t1), moving = oro2ants(readnii(muse)),
                                     transformlist = atlas_to_image$fwdtransforms, interpolator = "nearestNeighbor")
      muse_wmgm = antsApplyTransforms(fixed = oro2ants(t1), moving = oro2ants(readnii(wmgm)),
                                      transformlist = atlas_to_image$fwdtransforms, interpolator = "nearestNeighbor")
      antsImageWrite(muse_reg, paste0(output.dir_path_atlas, "/jlf_muse_reg",j, ".nii.gz"))
      antsImageWrite(muse_wmgm, paste0(output.dir_path_wmgm, "/jlf_muse_reg",j,"_wmgm.nii.gz"))
    }}
}
x= 1:length(dirs)
mclapply(x, do.jlf, mc.cores = 40) 
