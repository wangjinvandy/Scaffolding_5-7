# Scaffolding_5-7
 This is the code for paper "Neural representations of phonology in temporal cortex scaffold longitudinal reading gains in 5- to 7-year-old children" by Jin Wang, Marc F. Joanisse, James R. Booth
SPM12, ArtRepair, Marsbar toolbox are needed for these codes to run.

#Preprocessing:
main_phon.m is the main code for preprocessing. It calls the subfunctions in folder fmri_preproc_generic_mni, such as realignment, segmentation, normalization and smoothing

#Firstlevel analysis:
firstlevel_generate_phon.m is the first level analysis code. It calls the onsets, firstlevel_subfunctions folder. 

#Get betas:
make_paramObject.m and makeroi.m are two functions to make individual ROIs. You should run make_paramObject.m first to specify your parameters, and then type makeroi(paramObject) in the matlab command window to make the rois.
getbetas.m is the code to extract betas, by calling marsbar toolbox.

#Toolbox:
folder spm includes the functions we modified in our project. In spm_default.m, we changed defaults.mask.thresh =0.5 (instead of its default 0.8).
folder  art_repair includes the functions we used and modified in the preprocessing and firstlevel step to repair the volume with big movement and deweight these repaired images.