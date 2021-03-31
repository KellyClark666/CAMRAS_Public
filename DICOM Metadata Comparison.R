library(oro.nifti)
library(neurobase)
library(tidyverse)
library(dplyr)

#find differences between two csv's containg dicom metadata


Protocol_Fields = c("ProtocolName", "Field.Of.View", "Repetition.Time", "Echo.Time", 
                    "Flip.Angle", "Software.Versions", "Number.of.Averages", "Acquistion.Matrix", 
                    "Pixel.Spacing", "Slice.Thickness", "Scan.Time")                                #fields of interest

reference_hdr = read.csv("ref.csv")
target_hdr = read.csv("target.csv")

#fields of interest
Protocol_Fields = c("DICOM", "ProtocolName","SeriesInstanceUID", "Field.Of.View","Repetition.Time", "Echo.Time", 
                    "Flip.Angle", "Software.Versions", "Number.of.Averages", "Acquistion.Matrix", 
                    "Pixel.Spacing", "Slice.Thickness", "Scan.Time")                              

#subset ref DF based on specific protocol fields above
reference_hdr= reference_hdr[Protocol_Fields] 
reference_hdr = data.frame(lapply(reference_hdr, as.character), stringsAsFactors=FALSE)   #converts everything in DF to characters
reference_hdr[is.na(reference_hdr)]="Missing"  #replace NA with Missing, Missing = info not found in metadata

#find where PDT2_AX_TSE has an echo time of 11, and change protocol name accordingly
index_11 = which(reference_hdr$ProtocolName=="PDT2_AX_TSE" & reference_hdr$Echo.Time=="11" , arr.ind=TRUE)
for (i in 1:length(index_11)){
  reference_hdr[index_11[i],"ProtocolName"]="PDT2_AX_TSE_11"
}     

#find where PDT2_AX_TSE has an echo time of 101, and change protocol name accordingly
index_101 = which(reference_hdr$ProtocolName=="PDT2_AX_TSE" & reference_hdr$Echo.Time=="101" , arr.ind=TRUE)
for (i in 1:length(index_101)){
  reference_hdr[index_101[i],"ProtocolName"]="PDT2_AX_TSE_101"
}

#Above for loops are used to differentiate dual echo images PD and T2_AX_TSE since they have the same series instance ID

#get rid of duplicates
reference_hdr = reference_hdr[!duplicated(reference_hdr[c("SeriesInstanceUID","ProtocolName")]),]   
reference_hdr = reference_hdr %>% dplyr::rename_at((vars(-"ProtocolName")), paste0, "Reference_Patient") 
 
target_hdr = target_hdr[Protocol_Fields]
target_hdr = data.frame(lapply(target_hdr, as.character), stringsAsFactors=FALSE) 
target_hdr[is.na(target_hdr)]="Missing"
index_11 = which(target_hdr$ProtocolName=="PDT2_AX_TSE" & target_hdr$Echo.Time=="11" , arr.ind=TRUE)
for (i in 1:length(index_11)){
  target_hdr[index_11[i],"ProtocolName"]="PDT2_AX_TSE_11"
}

index_101 = which(target_hdr$ProtocolName=="PDT2_AX_TSE" & target_hdr$Echo.Time=="101" , arr.ind=TRUE)
for (i in 1:length(index_101)){
  target_hdr[index_101[i],"ProtocolName"]="PDT2_AX_TSE_101"
}

target_hdr = target_hdr[!duplicated(target_hdr[c("SeriesInstanceUID","ProtocolName")]),]
target_hdr = target_hdr %>% dplyr::rename_at((vars(-"ProtocolName")),paste0, "Trial_Patient")

#join and equate ref and target dataframes (equates by filling in missing data by duplicating existing data) 
merge = merge(reference_hdr, target_hdr, by ="ProtocolName")  

#subset merged DFs (ref and target) now that # of rows is equal between both
reference_Eq = merge %>% select("ProtocolName",ends_with("Reference_Patient"))   

#remove identifier added to column names so column names are equal between ref and target
names(reference_Eq) = gsub("Reference_Patient","",names(reference_Eq))      
reference_Eq = reference_Eq %>% arrange(ProtocolName,  Repetition.Time, Echo.Time,
                                  Flip.Angle, Software.Versions, Number.of.Averages, Acquistion.Matrix,
                                  Pixel.Spacing, Slice.Thickness)                              

target_Eq = merge %>% select("ProtocolName", ends_with("Trial_Patient"))
names(target_Eq) = gsub("Trial_Patient","",names(target_Eq))
target_Eq = target_Eq %>% arrange(ProtocolName ,Repetition.Time, Echo.Time,
                                    Flip.Angle, Software.Versions, Number.of.Averages, Acquistion.Matrix,
                                    Pixel.Spacing, Slice.Thickness)

Protocol_Fields_2 = c( "ProtocolName", "Field.Of.View", "Repetition.Time", "Echo.Time", 
                    "Flip.Angle", "Software.Versions", "Number.of.Averages", "Acquistion.Matrix", 
                    "Pixel.Spacing", "Slice.Thickness", "Scan.Time")

#subset to remove DICOM column for comparison since these will unequal between DFs and are irrelevant in comparison
sub_target_Eq = target_Eq[,-which(names(target_Eq)=="DICOM")]     
sub_target_Eq = sub_target_Eq[Protocol_Fields_2]

sub_ref_Eq=reference_Eq[,-which(names(reference_Eq)=="DICOM")]
sub_ref_Eq = sub_ref_Eq[Protocol_Fields_2]

#returns dataframe containing indices where differences can be found
diff = which(as.data.frame(sub_ref_Eq)!=as.data.frame(sub_target_Eq), arr.ind=TRUE)

row_index =as.numeric(diff[,1])       #save column containing row numbers as a list
column_index = as.numeric(diff[,2])   #save column containing column numbers as list
column = colnames(sub_target_Eq)      #saving column names to be indexed in for loop below

target_Eq= as.data.frame(target_Eq) %>% rename(DICOM_Patient=DICOM)       
reference_Eq = as.data.frame(reference_Eq) %>% rename(DICOM_Reference_Patient=DICOM)

data_DICOM = data.frame()    #create empty DF to be filled in by for loop
for(i in 1:length(row_index)){   #row_index (or col_index, the length will be the same) = number of differences between patient and healthy
  
  Modality= sub_target_Eq[row_index[i],1]                                  #Save modality name that corresponds to each difference 
  DICOM_Patient = as.character(target_Eq[row_index[i],2])                 #Save patient dicom file name where difference is coming from
  DICOM_Reference_Patient = as.character(reference_Eq[row_index[i],2])              
  Protocol.Field = column[column_index[i]]                                  #Save protocol field that corresponds to each difference
  Patient = sub_target_Eq[row_index[i], column_index[i]]                #find exact value of difference in TargetDiff based on index of difference
  Reference.Patient = sub_ref_Eq[row_index[i], column_index[i]] 
  bind = cbind(DICOM_Patient, DICOM_Reference_Patient,Modality,Protocol.Field, Patient, Reference.Patient) #column bind, modality, fields, patient value, and healthy value
  data_DICOM = rbind(data_DICOM, bind)                                                #row bind each iteration of bind to empty dataframe called data
  write.csv(data_DICOM, "Protocol_Differences.csv")               #write CSV of differences
  
}


