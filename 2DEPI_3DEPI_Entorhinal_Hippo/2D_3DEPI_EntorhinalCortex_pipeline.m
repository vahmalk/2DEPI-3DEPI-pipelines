%% Pipeline for assessing 2DEPI & 3DEPI data in Entorhinal cortex and Hippocampus
% - 2D EPI low res - high res
% - 3D EPI high res
% - tSNR, local TE and Jacobian
% - no cleaning vs cleaning via aCompCor
% - using Misun's pilot data 
% V Malekian, FIL Physics, 04/06/2024
clc
close all
clear all

addpath('/export/home/vmalekian/Malab_lib/NIfTI_20140122/')
addpath(genpath('/export/home/vmalekian/Yan/physio_correction/fmri_denoising-master'))


%% Structural preparations and segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir_inv2 = '/export/home/vmalekian/LC_fMRI/misun/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_largerFOV_INV2_0006/';
dir_struct = '/export/home/vmalekian/LC_fMRI/misun/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_largerFOV_UNI_Images_0008/';
unix(['fslreorient2std ' dir_struct 's*1.nii ' dir_struct 'T1_reo.nii.gz'])
unix(['fslreorient2std ' dir_inv2 's*1.nii ' dir_inv2 'T1_reo.nii.gz'])

unix(['fsl_anat -i ' dir_struct 'T1_reo -o ' dir_struct 'T1_reo --noreorient --nocrop'])
unix(['fsl_anat -i ' dir_inv2 'T1_reo -o ' dir_inv2 'T1_reo --noreorient --nocrop'])

% unix(['bet ' dir_inv2 'T1.anat/T1_biascorr ' dir_inv2 'T1.anat/T1_biascorr_bet  -f 0.5 -g 0'])
unix(['fslmaths ' dir_struct 'T1_reo.anat/T1_biascorr -mul ' dir_inv2 'T1_reo.anat/T1_biascorr_brain_mask ' dir_struct 'T1_reo.anat/T1_biascorr_brain_final'])


%SPM tissue segmentation
unix(['gunzip ' (dir_struct) 'T1_reo.anat/T1_biascorr.nii.gz'])

matlabbatch{1}.spm.spatial.preproc.channel.vols = {
    [dir_struct 'T1_reo.anat/T1_biascorr.nii,1']
    }

matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.01;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 30;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 2;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.0005 0.25 0.025 0.1];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 1;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
    NaN NaN NaN];

spm_jobman('run', matlabbatch);

unix(['fslmaths ' dir_struct 'T1_reo.anat/c1T1_biascorr -thr 0.5 -bin ' dir_struct 'c1T1_reo_bin']) %% WM
unix(['fslmaths ' dir_struct 'T1_reo.anat/c2T1_biascorr -thr 0.5 -bin ' dir_struct 'c2T1_reo_bin']) %% WM
unix(['fslmaths ' dir_struct 'T1_reo.anat/c3T1_biascorr -thr 0.5 -bin ' dir_struct 'c3T1_reo_bin']) %% WM

unix(['fslmaths ' dir_struct 'c1T1_reo_bin -add ' dir_struct 'c2T1_reo_bin -add ' dir_struct 'c3T1_reo_bin ' dir_struct 'c_reo_bin'])
unix(['fslmaths ' dir_struct 'c_reo_bin -div 3 ' dir_struct 'c_reo_bin_mask'])
unix(['fslmaths ' dir_struct 'c_reo_bin_mask -fillh ' dir_struct 'c_reo_bin_mask_fill'])

unix(['fslmaths ' dir_struct 'c_reo_bin_mask -dilF ' dir_struct 'c_reo_bin_mask_dilF'])
unix(['fslmaths ' dir_struct 'c_reo_bin_mask_dilF -eroF ' dir_struct 'c_reo_bin_mask_dilF_eroF'])
unix(['fslmaths ' dir_struct 'c_reo_bin_mask_dilF_eroF -fillh ' dir_struct 'c_reo_bin_mask_dilF_eroF_fill'])
unix(['fslmaths ' dir_struct 'T1_reo -mul ' dir_struct 'c_reo_bin_mask_dilF_eroF_fill ' dir_struct 'T1_reo_brain'])
% unix(['fsl_anat -i ' dir_struct 'T1_reo_brain -o ' dir_struct 'T1_reo_brain --noreorient --nocrop'])


%% masking
unix(['bet ' dir_epi1_feat 'mean_func ' dir_epi1_feat 'mean_func_brain -f 0.3 -g 0 -m'])
% unix(['bet ' dir_epi2_feat 'mean_func ' dir_epi2_feat 'mean_func_brain -f 0.075 -g 0 -m'])

unix(['fslmaths ' dir_epi1_feat 'filtered_func_data_snr.nii.gz -mas ' dir_epi1_feat 'mean_func_brain_mask ' dir_epi1_feat 'filtered_func_data_snr_masked'])
unix(['fslmaths ' dir_epi1_feat 'filtered_func_data_cleaned_snr.nii.gz -mas ' dir_epi1_feat 'mean_func_brain_mask ' dir_epi1_feat 'filtered_func_data_cleaned_snr_masked'])

%% Start of the pipeline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dir_epi1 = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/cmrr_mbep2d_1point25_hiedi_GRE_108_AP_0014/';
% dir_epi2 = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/cmrr_mbep2d_1point25_hiedi_GRE_108_PA_0016/';
% dir_epi1_ref = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/cmrr_mbep2d_1point25_hiedi_GRE_108_AP_SBRef_0013/';
% dir_epi2_ref = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/cmrr_mbep2d_1point25_hiedi_GRE_108_PA_SBRef_0015/';

dir_epi1 = '/export/home/vmalekian/LC_fMRI/misun/cmrr_mbep2d_1p25_fMRI_PA_0012/';
dir_epi2 = '/export/home/vmalekian/LC_fMRI/misun/cmrr_mbep2d_1p00_fMRI_PA_0016/';
dir_epi3 = '/export/home/vmalekian/LC_fMRI/misun/nc_epi3d_0p92mm_noSEG_BWT18_180_300us_0021/';


dir_epi1_rev = '/export/home/vmalekian/LC_fMRI/misun/cmrr_mbep2d_1p25_REV_0010/';
dir_epi2_rev = '/export/home/vmalekian/LC_fMRI/misun/cmrr_mbep2d_1p00_REV_0014/';

dir_distortion_1 = '/export/home/vmalekian/LC_fMRI/misun/distortion_topup_12mm/';
dir_distortion_2 = '/export/home/vmalekian/LC_fMRI/misun/distortion_topup_10mm/';
dir_distortion_3 = '/export/home/vmalekian/LC_fMRI/misun/distortion_topup_09mm/';


unix(['fslmerge -t ' dir_epi1 'func ' dir_epi1 'f*.nii'])
unix(['fslmerge -tr ' dir_epi1 'func ' dir_epi1 'func 2.1'])

unix(['fslmerge -t ' dir_epi2 'func ' dir_epi2 'f*.nii'])
unix(['fslmerge -tr ' dir_epi2 'func ' dir_epi2 'func 2.6'])

unix(['fslmerge -t ' dir_epi3 'func ' dir_epi3 'f*.nii'])
unix(['fslmerge -tr ' dir_epi3 'func ' dir_epi3 'func 2.6'])

%%MC and merging PA and AP for 0.92 mm, 4 vol
unix(['fslroi ' dir_epi3 'func ' dir_epi3 'AP_fst 0 -1 0 -1 0 -1 0 4'])
unix(['fslroi ' dir_epi3 'func ' dir_epi3 'func_main 0 -1 0 -1 0 -1 4 122'])

unix(['mcflirt -in ' dir_epi3 'AP_fst -o ' dir_epi3 'AP_fst_mc -refvol 3 -spline_final -plots -report -stats -mats'])
unix(['fslmaths ' dir_epi3 'AP_fst_mc -thr 0 ' dir_epi3 'AP_fst_mc'])
unix(['fslmaths ' dir_epi3 'AP_fst_mc -Tmean ' dir_epi3 'AP_tmean'])

unix(['mcflirt -in ' dir_epi3 'func_main -o ' dir_epi3 'func_main_mc -refvol 0 -spline_final -plots -report -stats -mats'])
unix(['fslmaths ' dir_epi3 'func_main_mc -thr 0 ' dir_epi3 'func_main_mc'])
unix(['fslroi ' dir_epi3 'func_main_mc ' dir_epi3 'PA_fst 0 -1 0 -1 0 -1 0 4'])
unix(['fslmaths ' dir_epi3 'PA_fst -Tmean ' dir_epi3 'PA_tmean'])

unix(['fslmerge -t ' dir_distortion_3 'PA_AP ' dir_epi3 'PA_tmean ' dir_epi3 'AP_tmean'])


%%MC and merging PA and AP for 1.25 mm, 1 vol
unix(['fslmerge -t ' dir_epi1_rev 'func ' dir_epi1_rev 'f*.nii'])
unix(['fslmerge -tr ' dir_epi1_rev 'func ' dir_epi1_rev 'func 2.1'])
unix(['mcflirt -in ' dir_epi1_rev 'func -o ' dir_epi1_rev 'func_mc -refvol 3 -spline_final -plots -report -stats -mats'])
unix(['fslmaths ' dir_epi1_rev 'func_mc -thr 0 ' dir_epi1_rev 'func_mc'])
unix(['fslmaths ' dir_epi1_rev 'func_mc -Tmean ' dir_epi1_rev 'AP_tmean'])

unix(['mcflirt -in ' dir_epi1 'func -o ' dir_epi1 'func_mc -refvol 0 -spline_final -plots -report -stats -mats'])
unix(['fslmaths ' dir_epi1 'func_mc -thr 0 ' dir_epi1 'func_mc'])
unix(['fslroi ' dir_epi1 'func_mc ' dir_epi1 'PA_fst 0 -1 0 -1 0 -1 0 4'])
unix(['fslmaths ' dir_epi1 'PA_fst -Tmean ' dir_epi1 'PA_tmean'])

unix(['fslmerge -t ' dir_distortion_1 'PA_AP ' dir_epi1 'PA_tmean ' dir_epi1_rev 'AP_tmean'])

%%MC and merging PA and AP for 1 mm, 4 vols
unix(['fslmerge -t ' dir_epi2_rev 'func ' dir_epi2_rev 'f*.nii'])
unix(['fslmerge -tr ' dir_epi2_rev 'func ' dir_epi2_rev 'func 2.6'])
unix(['mcflirt -in ' dir_epi2_rev 'func -o ' dir_epi2_rev 'func_mc -refvol 3 -spline_final -plots -report -stats -mats'])
unix(['fslmaths ' dir_epi2_rev 'func_mc -thr 0 ' dir_epi2_rev 'func_mc'])
unix(['fslmaths ' dir_epi2_rev 'func_mc -Tmean ' dir_epi2_rev 'AP_tmean'])

unix(['mcflirt -in ' dir_epi2 'func -o ' dir_epi2 'func_mc -refvol 0 -spline_final -plots -report -stats -mats'])
unix(['fslmaths ' dir_epi2 'func_mc -thr 0 ' dir_epi2 'func_mc'])
unix(['fslroi ' dir_epi2 'func_mc ' dir_epi2 'PA_fst 0 -1 0 -1 0 -1 0 4'])
unix(['fslmaths ' dir_epi2 'PA_fst -Tmean ' dir_epi2 'PA_tmean'])

unix(['fslmerge -t ' dir_distortion_2 'PA_AP ' dir_epi2 'PA_tmean ' dir_epi2_rev 'AP_tmean'])


%%Effective Echo spacing
%%1.00 mm 2D-EPI= 780*197/3 = 51.220 ms
%%1.25 mm 2D-EPI= 780*159/3 = 41.340 ms
%%0.92 mm 3D-EPI= 1220*207/4 = 63.135 ms

%%topup 1.25mm 2DEPI
unix(['topup --imain=' dir_distortion_1 'PA_AP --datain=' dir_distortion_1 'param.txt --config=b02b0.cnf --out=' dir_distortion_1 'topup_PA_AP --fout=' dir_distortion_1 'topup_field ' ...
    '--iout=' dir_distortion_1 'topup_unwarp --jacout=' dir_distortion_1 'jac --rbmout=' dir_distortion_1 'xfm --dfout=' dir_distortion_1 'topup_warpfield'])
unix(['applytopup --imain=' dir_epi1 'func_mc --topup=' dir_distortion_1 'topup_PA_AP --datain=' dir_distortion_1 'param.txt --inindex=1 --out=' dir_epi1 'func_mc_dc --method=jac']) %% applytopup with Jacobian correction
unix(['fslmaths ' dir_epi1 'func_mc_dc -thr 0 ' dir_epi1 'func_mc_dc'])

%%topup 1.00mm 2DEPI
unix(['topup --imain=' dir_distortion_2 'PA_AP --datain=' dir_distortion_2 'param.txt --config=b02b0.cnf --out=' dir_distortion_2 'topup_PA_AP --fout=' dir_distortion_2 'topup_field ' ...
    '--iout=' dir_distortion_2 'topup_unwarp --jacout=' dir_distortion_2 'jac --rbmout=' dir_distortion_2 'xfm --dfout=' dir_distortion_2 'topup_warpfield'])
unix(['applytopup --imain=' dir_epi2 'func_mc --topup=' dir_distortion_2 'topup_PA_AP --datain=' dir_distortion_2 'param.txt --inindex=1 --out=' dir_epi2 'func_mc_dc --method=jac']) %% applytopup with Jacobian correction
unix(['fslmaths ' dir_epi2 'func_mc_dc -thr 0 ' dir_epi2 'func_mc_dc'])

%%topup 0.92mm 3DEPI
unix(['topup --imain=' dir_distortion_3 'PA_AP --datain=' dir_distortion_3 'param.txt --config=b02b0.cnf --out=' dir_distortion_3 'topup_PA_AP --fout=' dir_distortion_3 'topup_field ' ...
    '--iout=' dir_distortion_3 'topup_unwarp --jacout=' dir_distortion_3 'jac --rbmout=' dir_distortion_3 'xfm --dfout=' dir_distortion_3 'topup_warpfield'])
unix(['applytopup --imain=' dir_epi3 'func_main_mc --topup=' dir_distortion_3 'topup_PA_AP --datain=' dir_distortion_3 'param.txt --inindex=1 --out=' dir_epi3 'func_mc_dc --method=jac']) %% applytopup with Jacobian correction
unix(['fslmaths ' dir_epi3 'func_mc_dc -thr 0 ' dir_epi3 'func_mc_dc'])

%%topup 0.92mm for WB 3DEPI
dir_epi3_wb = '/export/home/vmalekian/LC_fMRI/misun/nc_epi3d_0p92mm_noSEG_WB_0024/';
dir_distortion_3_wb = '/export/home/vmalekian/LC_fMRI/misun/distortion_topup_09mm_wb/';

unix(['fslmerge -t ' dir_distortion_3_wb 'PA_AP ' dir_epi3_wb 'PA ' dir_epi3_wb 'AP'])

unix(['topup --imain=' dir_distortion_3_wb 'PA_AP --datain=' dir_distortion_3_wb 'param.txt --config=b02b0.cnf --out=' dir_distortion_3_wb 'topup_PA_AP --fout=' dir_distortion_3_wb 'topup_field ' ...
    '--iout=' dir_distortion_3_wb 'topup_unwarp --jacout=' dir_distortion_3_wb 'jac --rbmout=' dir_distortion_3_wb 'xfm --dfout=' dir_distortion_3_wb 'topup_warpfield'])
unix(['applytopup --imain=' dir_epi3_wb 'PA --topup=' dir_distortion_3_wb 'topup_PA_AP --datain=' dir_distortion_3_wb 'param.txt --inindex=1 --out=' dir_epi3_wb 'PA_dc --method=jac']) %% applytopup with Jacobian correction
unix(['fslmaths ' dir_epi3_wb 'PA_dc -thr 0 ' dir_epi3_wb 'PA_dc'])

%% Perform HP, BET, Co-registeration, Normalisation
unix(['fslmaths ' dir_epi1 'func_mc_dc -Tmean ' dir_epi1 'func_mc_dc_tmean'])
unix(['fslmaths ' dir_epi2 'func_mc_dc -Tmean ' dir_epi2 'func_mc_dc_tmean'])
unix(['fslmaths ' dir_epi3 'func_mc_dc -Tmean ' dir_epi3 'func_mc_dc_tmean'])

unix (['fslmaths ' dir_epi1 'func_mc_dc -bptf 23.809 -1 -add ' dir_epi1 'func_mc_dc_tmean ' dir_epi1 'func_mc_dc_hp']) % 1/(2*0.01*2.1)
unix (['fslmaths ' dir_epi2 'func_mc_dc -bptf 19.231 -1 -add ' dir_epi2 'func_mc_dc_tmean ' dir_epi2 'func_mc_dc_hp']) % 1/(2*0.01*2.6)
unix (['fslmaths ' dir_epi3 'func_mc_dc -bptf 19.231 -1 -add ' dir_epi3 'func_mc_dc_tmean ' dir_epi3 'func_mc_dc_hp']) % 1/(2*0.01*2.6)

unix(['fslmaths ' dir_epi1 'func_mc_dc_hp -thr 0 ' dir_epi1 'func_mc_dc_hp'])
unix(['fslmaths ' dir_epi2 'func_mc_dc_hp -thr 0 ' dir_epi2 'func_mc_dc_hp'])
unix(['fslmaths ' dir_epi3 'func_mc_dc_hp -thr 0 ' dir_epi3 'func_mc_dc_hp'])

%%3DEPI slice removal
unix(['fslroi ' dir_epi3 'func_mc_dc_hp ' dir_epi3 'func_mc_dc_hpc 0 -1 0 -1 5 78 0 -1'])
unix(['fslroi ' dir_epi3 'func_mc_dc_tmean ' dir_epi3 'func_mc_dc_tmeanc 0 -1 0 -1 5 78'])

%% Structral registeration

unix(['epi_reg --epi=' dir_epi1 'func_mc_dc_tmean --t1=' dir_struct 'T1_reo.anat/T1_biascorr --t1brain=' dir_struct 'T1_reo.anat/T1_biascorr_brain --out=' dir_epi1 'func_mc_dc_tmean_highres --pedir=y'])
unix(['epi_reg --epi=' dir_epi2 'func_mc_dc_tmean --t1=' dir_struct 'T1_reo.anat/T1_biascorr --t1brain=' dir_struct 'T1_reo.anat/T1_biascorr_brain --out=' dir_epi2 'func_mc_dc_tmean_highres --pedir=y'])

unix(['epi_reg --epi=' dir_epi3_wb 'PA_dc --t1=' dir_struct 'T1_reo.anat/T1_biascorr --t1brain=' dir_struct 'T1_reo.anat/T1_biascorr_brain --out=' dir_epi3_wb 'PA_dc_highres --pedir=y'])
unix(['flirt -in ' dir_epi3 'func_mc_dc_tmeanc -ref ' dir_epi3_wb 'PA_dc -out ' dir_epi3 'func_mc_dc_tmean_wb -omat ' dir_epi3 'func2wb.mat -dof 6 -interp spline'])
unix(['convert_xfm -omat ' dir_epi3 'func_mc_dc_tmean_highres.mat -concat ' dir_epi3_wb 'PA_dc_highres.mat ' dir_epi3 'func2wb.mat'])
unix(['applywarp --ref=' dir_struct 'T1_reo.anat/T1_biascorr_brain --in=' dir_epi3 'func_mc_dc_tmeanc --out=' dir_epi3 'func_mc_dc_tmean_highres --premat=' dir_epi3 'func_mc_dc_tmean_highres.mat --interp=spline'])

%% aCompCor CSF and WM
unix(['fslmaths ' dir_struct 'T1_reo.anat/c2T1_biascorr -thr 0.998 -bin ' dir_struct 'c2T1_reo_mask']) %% WM
unix(['fslmaths ' dir_struct 'T1_reo.anat/c3T1_biascorr -thr 0.998 -bin ' dir_struct 'c3T1_reo_mask']) %% CSF

unix(['bet ' dir_epi1 'func_mc_dc_tmean ' dir_epi1 'func_mc_dc_tmean_brain -f 0.5 -g 0 -m'])
unix(['bet ' dir_epi2 'func_mc_dc_tmean ' dir_epi2 'func_mc_dc_tmean_brain -f 0.5 -g 0 -m'])
unix(['bet ' dir_epi3 'func_mc_dc_tmeanc ' dir_epi3 'func_mc_dc_tmean_brain -f 0.5 -g 0 -m'])

unix(['convert_xfm -inverse -omat ' dir_epi1 'func_mc_dc_tmean_highres_inv.mat ' dir_epi1 'func_mc_dc_tmean_highres.mat'])
unix(['convert_xfm -inverse -omat ' dir_epi2 'func_mc_dc_tmean_highres_inv.mat ' dir_epi2 'func_mc_dc_tmean_highres.mat'])
unix(['convert_xfm -inverse -omat ' dir_epi3 'func_mc_dc_tmean_highres_inv.mat ' dir_epi3 'func_mc_dc_tmean_highres.mat'])


unix(['applywarp --ref=' dir_epi1 'func_mc_dc_tmean --in=' dir_struct 'c2T1_reo_mask' ...
    ' --out=' dir_epi1 'c2T1_reo_2fMREPI --premat=' dir_epi1 'func_mc_dc_tmean_highres_inv.mat --interp=nn'])
unix(['applywarp --ref=' dir_epi1 'func_mc_dc_tmean --in=' dir_struct 'c3T1_reo_mask' ...
    ' --out=' dir_epi1 'c3T1_reo_2fMREPI --premat=' dir_epi1 'func_mc_dc_tmean_highres_inv.mat --interp=nn'])

unix(['applywarp --ref=' dir_epi2 'func_mc_dc_tmean --in=' dir_struct 'c2T1_reo_mask' ...
    ' --out=' dir_epi2 'c2T1_reo_2fMREPI --premat=' dir_epi2 'func_mc_dc_tmean_highres_inv.mat --interp=nn'])
unix(['applywarp --ref=' dir_epi2 'func_mc_dc_tmean --in=' dir_struct 'c3T1_reo_mask' ...
    ' --out=' dir_epi2 'c3T1_reo_2fMREPI --premat=' dir_epi2 'func_mc_dc_tmean_highres_inv.mat --interp=nn'])

unix(['applywarp --ref=' dir_epi3 'func_mc_dc_tmeanc --in=' dir_struct 'c2T1_reo_mask' ...
    ' --out=' dir_epi3 'c2T1_reo_2fMREPI --premat=' dir_epi3 'func_mc_dc_tmean_highres_inv.mat --interp=nn'])
unix(['applywarp --ref=' dir_epi3 'func_mc_dc_tmeanc --in=' dir_struct 'c3T1_reo_mask' ...
    ' --out=' dir_epi3 'c3T1_reo_2fMREPI --premat=' dir_epi3 'func_mc_dc_tmean_highres_inv.mat --interp=nn'])


unix(['fslmaths ' dir_epi1 'c2T1_reo_2fMREPI -mul ' dir_epi1 'func_mc_dc_tmean_brain_mask ' dir_epi1 'c2T1_reo_2fMREPI_BET'])
unix(['fslmaths ' dir_epi1 'c3T1_reo_2fMREPI -mul ' dir_epi1 'func_mc_dc_tmean_brain_mask ' dir_epi1 'c3T1_reo_2fMREPI_BET'])

unix(['fslmaths ' dir_epi2 'c2T1_reo_2fMREPI -mul ' dir_epi2 'func_mc_dc_tmean_brain_mask ' dir_epi2 'c2T1_reo_2fMREPI_BET'])
unix(['fslmaths ' dir_epi2 'c3T1_reo_2fMREPI -mul ' dir_epi2 'func_mc_dc_tmean_brain_mask ' dir_epi2 'c3T1_reo_2fMREPI_BET'])

unix(['fslmaths ' dir_epi3 'c2T1_reo_2fMREPI -mul ' dir_epi3 'func_mc_dc_tmean_brain_mask ' dir_epi3 'c2T1_reo_2fMREPI_BET'])
unix(['fslmaths ' dir_epi3 'c3T1_reo_2fMREPI -mul ' dir_epi3 'func_mc_dc_tmean_brain_mask ' dir_epi3 'c3T1_reo_2fMREPI_BET'])


%%2DEPI 1.25mm physio regressors
string = 'c2T1_reo_2fMREPI_BET.nii.gz';
a1m=dir([dir_epi1,string]);
nf=load_untouch_nii([dir_epi1,a1m(1).name]);
roi_W= (nf.img);
string = 'c3T1_reo_2fMREPI_BET.nii.gz';
a1m=dir([dir_epi1,string]);
nf=load_untouch_nii([dir_epi1,a1m(1).name]);
roi_C= (nf.img);

string = 'func_mc_dc_hp.nii.gz';
a1m=dir([dir_epi1,string]);
nf=load_untouch_nii([dir_epi1,a1m(1).name]);
data= (nf.img);
rois1 = {roi_C,roi_W};
X= fmri_compcor(data,rois1,[1/4 1/4]);% Alternitvely can use X= fmri_compcor(data,rois,[0.25 0.25])

TR=2.1;
t = 0:TR:315; % time vector
Fs = 1/(t(2)-t(1)); % Sampling frequency (Hz)
PxxT=zeros(floor(size(X,1)/2),size(X,2));

%%remove any regressers with low frequency contents using power spectrum feature
for fk=1:size(X,2)
    x = X(:,fk);
    N = length(x); % Length of the signal
    Xf = fft(x,N); % Compute FFT
    Pxx = abs(Xf).^2/N; % Power spectrum
    f = (0:N-1)*(Fs/N); % Frequency vector
    Pxx1=Pxx(1:length(Pxx)/2);
    f1=f(1:length(f)/2);
    PxxT(:,fk) = Pxx1;
end
[np,mp] = size(PxxT);

freqfeat = sum(PxxT(1:round(np/2),:))./sum(PxxT(1:round(np/1),:));
final_X = X(:,find(freqfeat<(1/2)==1));


for i2= 1:size(final_X, 2)
    filename = [dir_epi1 'physio_data_' sprintf('%02d',i2) '.txt'];
    fileID = fopen(filename,'w');

    for i1 = 1:size(final_X, 1)
        fprintf(fileID, '%f ', final_X(i1, i2));
        fprintf(fileID, '\n');
    end

    fclose(fileID);
end
writematrix(final_X,[dir_epi1 'physio.txt'])


[res,model_info] = fmri_cleaning(data,-1,final_X);
string_compcor = [dir_epi1 'func_mc_dc_hp_physio.nii.gz'];
save_nii(make_nii((res)),string_compcor)
unix(['fslcpgeom ' dir_epi1 'func_mc_dc_hp.nii.gz ' string_compcor])

%%2DEPI 1mm physio regressors
string = 'c2T1_reo_2fMREPI_BET.nii.gz';
a1m=dir([dir_epi2,string]);
nf=load_untouch_nii([dir_epi2,a1m(1).name]);
roi_W= (nf.img);
string = 'c3T1_reo_2fMREPI_BET.nii.gz';
a1m=dir([dir_epi2,string]);
nf=load_untouch_nii([dir_epi2,a1m(1).name]);
roi_C= (nf.img);
string = 'func_mc_dc_hp.nii.gz';
a1m=dir([dir_epi2,string]);
nf=load_untouch_nii([dir_epi2,a1m(1).name]);
data= (nf.img);
rois2 = {roi_C,roi_W};
X= fmri_compcor(data,rois2,[1/4 1/4]);% Alternitvely can use X= fmri_compcor(data,rois,[0.25 0.25])

%%remove any regressers with low frequency contents using power spectrum feature
TR=2.6;
t = 0:TR:315; % time vector
Fs = 1/(t(2)-t(1)); % Sampling frequency (Hz)
PxxT=zeros(floor(size(X,1)/2),size(X,2));

for fk=1:size(X,2)
    x = X(:,fk);
    N = length(x); % Length of the signal
    Xf = fft(x,N); % Compute FFT
    Pxx = abs(Xf).^2/N; % Power spectrum
    f = (0:N-1)*(Fs/N); % Frequency vector
    Pxx1=Pxx(1:length(Pxx)/2);
    f1=f(1:length(f)/2);
    PxxT(:,fk) = Pxx1;
end
[np,mp] = size(PxxT);

freqfeat = sum(PxxT(1:round(np/2),:))./sum(PxxT(1:round(np/1),:));
final_X = X(:,find(freqfeat<(1/2)==1));

for i2= 1:size(final_X, 2)
    filename = [dir_epi2 'physio_data_' sprintf('%02d',i2) '.txt'];
    fileID = fopen(filename,'w');

    for i1 = 1:size(final_X, 1)
        fprintf(fileID, '%f ', final_X(i1, i2));
        fprintf(fileID, '\n');
    end

    fclose(fileID);
end

writematrix(final_X,[dir_epi2 'physio.txt'])


[res,model_info] = fmri_cleaning(data,-1,final_X);
string_compcor = [dir_epi2 'func_mc_dc_hp_physio.nii.gz'];
save_nii(make_nii((res)),string_compcor)
unix(['fslcpgeom ' dir_epi2 'func_mc_dc_hp.nii.gz ' string_compcor])

%%3DEPI 0.9mm physio regressors
string = 'c2T1_reo_2fMREPI_BET.nii.gz';
a1m=dir([dir_epi3,string]);
nf=load_untouch_nii([dir_epi3,a1m(1).name]);
roi_W= (nf.img);
string = 'c3T1_reo_2fMREPI_BET.nii.gz';
a1m=dir([dir_epi3,string]);
nf=load_untouch_nii([dir_epi3,a1m(1).name]);
roi_C= (nf.img);
string = 'func_mc_dc_hpc.nii.gz';
a1m=dir([dir_epi3,string]);
nf=load_untouch_nii([dir_epi3,a1m(1).name]);
data= (nf.img);
rois2 = {roi_C,roi_W};
X= fmri_compcor(data,rois2,[1/4 1/4]);% Alternitvely can use X= fmri_compcor(data,rois,[0.25 0.25])

%%remove any regressers with low frequency contents using power spectrum feature
TR=2.6;
t = 0:TR:315; % time vector
Fs = 1/(t(2)-t(1)); % Sampling frequency (Hz)
PxxT=zeros(floor(size(X,1)/2),size(X,2));

for fk=1:size(X,2)
    x = X(:,fk);
    N = length(x); % Length of the signal
    Xf = fft(x,N); % Compute FFT
    Pxx = abs(Xf).^2/N; % Power spectrum
    f = (0:N-1)*(Fs/N); % Frequency vector
    Pxx1=Pxx(1:length(Pxx)/2);
    f1=f(1:length(f)/2);
    PxxT(:,fk) = Pxx1;
end
[np,mp] = size(PxxT);

freqfeat = sum(PxxT(1:round(np/2),:))./sum(PxxT(1:round(np/1),:));
final_X = X(:,find(freqfeat<(1/2)==1));


for i2= 1:size(final_X, 2)
    filename = [dir_epi3 'physio_data_' sprintf('%02d',i2) '.txt'];
    fileID = fopen(filename,'w');

    for i1 = 1:size(final_X, 1)
        fprintf(fileID, '%f ', final_X(i1, i2));
        fprintf(fileID, '\n');
    end

    fclose(fileID);
end

writematrix(final_X,[dir_epi3 'physio.txt'])

[res,model_info] = fmri_cleaning(data,-1,final_X);
string_compcor = [dir_epi3 'func_mc_dc_hpc_physio.nii.gz'];
save_nii(make_nii((res)),string_compcor)
unix(['fslcpgeom ' dir_epi3 'func_mc_dc_hpc.nii.gz ' string_compcor])
 

%% tsnr calculations %%

unix(['fslmaths ' dir_epi1 'func_mc_dc_hp -Tmean ' dir_epi1 'func_mc_dc_hp_tmean'])
unix(['fslmaths ' dir_epi1 'func_mc_dc_hp -Tstd ' dir_epi1 'func_mc_dc_hp_std'])
unix(['fslmaths ' dir_epi1 'func_mc_dc_hp_tmean -div ' dir_epi1 'func_mc_dc_hp_std ' dir_epi1 'func_mc_dc_hp_snr'])

unix(['fslmaths ' dir_epi2 'func_mc_dc_hp -Tmean ' dir_epi2 'func_mc_dc_hp_tmean'])
unix(['fslmaths ' dir_epi2 'func_mc_dc_hp -Tstd ' dir_epi2 'func_mc_dc_hp_std'])
unix(['fslmaths ' dir_epi2 'func_mc_dc_hp_tmean -div ' dir_epi2 'func_mc_dc_hp_std ' dir_epi2 'func_mc_dc_hp_snr'])

unix(['fslmaths ' dir_epi3 'func_mc_dc_hpc -Tmean ' dir_epi3 'func_mc_dc_hp_tmean'])
unix(['fslmaths ' dir_epi3 'func_mc_dc_hpc -Tstd ' dir_epi3 'func_mc_dc_hp_std'])
unix(['fslmaths ' dir_epi3 'func_mc_dc_hp_tmean -div ' dir_epi3 'func_mc_dc_hp_std ' dir_epi3 'func_mc_dc_hp_snr'])


%% tsnr calculations after cleaning with aCompCor
unix(['fslmaths ' dir_epi1 'func_mc_dc_hp_physio -Tmean ' dir_epi1 'func_mc_dc_hp_physio_tmean'])
unix(['fslmaths ' dir_epi1 'func_mc_dc_hp_physio -Tstd ' dir_epi1 'func_mc_dc_hp_physio_std'])
unix(['fslmaths ' dir_epi1 'func_mc_dc_hp_physio_tmean -div ' dir_epi1 'func_mc_dc_hp_physio_std ' dir_epi1 'func_mc_dc_hp_physio_snr'])

unix(['fslmaths ' dir_epi2 'func_mc_dc_hp_physio -Tmean ' dir_epi2 'func_mc_dc_hp_physio_tmean'])
unix(['fslmaths ' dir_epi2 'func_mc_dc_hp_physio -Tstd ' dir_epi2 'func_mc_dc_hp_physio_std'])
unix(['fslmaths ' dir_epi2 'func_mc_dc_hp_physio_tmean -div ' dir_epi2 'func_mc_dc_hp_physio_std ' dir_epi2 'func_mc_dc_hp_physio_snr'])

unix(['fslmaths ' dir_epi3 'func_mc_dc_hpc_physio -Tmean ' dir_epi3 'func_mc_dc_hpc_physio_tmean'])
unix(['fslmaths ' dir_epi3 'func_mc_dc_hpc_physio -Tstd ' dir_epi3 'func_mc_dc_hpc_physio_std'])
unix(['fslmaths ' dir_epi3 'func_mc_dc_hpc_physio_tmean -div ' dir_epi3 'func_mc_dc_hpc_physio_std ' dir_epi3 'func_mc_dc_hpc_physio_snr'])


%% MNI_2009C atlas normalisation for extracting brainstem

darya_atlas_09c = '/export/home/vmalekian/LC_fMRI/20231207.M700809_FIL/darya_atlas/mni_icbm152_nlin_asym_09c/';
unix(['fslmaths ' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c -mas ' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c_mask ' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c_masked'])

unix(['flirt -in ' dir_struct 'T1_reo.anat/T1_biascorr_brain -ref ' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c_masked -out ' dir_struct 'highres2mni9c -omat ' dir_struct 'highres2mni9c.mat -cost corratio -dof 12 -interp spline'])
unix(['fnirt --in=' dir_struct 'T1_reo.anat/T1_biascorr --ref=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c --aff=' dir_struct 'highres2mni9c.mat' ...
    ' --cout=' dir_struct 'highres2mni9c_warp --iout=' dir_struct 'highres2mni9c_out --fout=' dir_struct 'highres2mni9c_field --interp=spline']) 


unix(['flirt -in ' dir_struct 'T1_reo.anat/T1_biascorr_brain -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain -out ' dir_struct 'highres2mni1mm -omat ' dir_struct 'highres2mni1mm.mat -cost corratio -dof 12 -interp spline'])
unix(['fnirt --in=' dir_struct 'T1_reo.anat/T1_biascorr --ref=/usr/local/fsl/data/standard/MNI152_T1_1mm --aff=' dir_struct 'highres2mni1mm.mat' ...
    ' --cout=' dir_struct 'highres2mni1mm_warp --iout=' dir_struct 'highres2mni1mm_out --fout=' dir_struct 'highres2mni1mm_field --interp=spline']) 

%% Create transformation matrices
unix(['convertwarp --ref=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c --premat=' dir_epi1 'func_mc_dc_tmean_highres.mat --warp1=' dir_struct 'highres2mni9c_warp --out=' dir_epi1 'func2atlas_warp'])
unix(['invwarp -w ' dir_epi1 'func2atlas_warp -o ' dir_epi1 'func2atlas_warp_inv -r ' dir_epi1 'func_mc_dc_tmean'])
unix(['applywarp --ref=' dir_epi1 'func_mc_dc_tmean --in=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c_masked --out=' dir_epi1 'mni_icbm152_t1_tal_nlin_asym_09c_masked2fmri --warp=' dir_epi1 'func2atlas_warp_inv --interp=spline'])

unix(['convertwarp --ref=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c --premat=' dir_epi2 'func_mc_dc_tmean_highres.mat --warp1=' dir_struct 'highres2mni9c_warp --out=' dir_epi2 'func2atlas_warp'])
unix(['invwarp -w ' dir_epi2 'func2atlas_warp -o ' dir_epi2 'func2atlas_warp_inv -r ' dir_epi2 'func_mc_dc_tmean'])
unix(['applywarp --ref=' dir_epi2 'func_mc_dc_tmean --in=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c_masked --out=' dir_epi2 'mni_icbm152_t1_tal_nlin_asym_09c_masked2fmri --warp=' dir_epi2 'func2atlas_warp_inv --interp=spline'])

unix(['convertwarp --ref=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c --premat=' dir_epi3 'func_mc_dc_tmean_highres.mat --warp1=' dir_struct 'highres2mni9c_warp --out=' dir_epi3 'func2atlas_warp'])
unix(['invwarp -w ' dir_epi3 'func2atlas_warp -o ' dir_epi3 'func2atlas_warp_inv -r ' dir_epi3 'func_mc_dc_tmeanc'])
unix(['applywarp --ref=' dir_epi3 'func_mc_dc_tmeanc --in=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c_masked --out=' dir_epi3 'mni_icbm152_t1_tal_nlin_asym_09c_masked2fmri --warp=' dir_epi3 'func2atlas_warp_inv --interp=spline'])


%% MNI_2009C Hippocampus mask labels 
rhp = 48;
lhp = 99;
string_ats= 'mni_icbm152_CerebrA_tal_nlin_sym_09c.nii'; 
a0m=dir([darya_atlas_09c,string_ats]);
nf=load_untouch_nii([darya_atlas_09c,a0m(1).name]);
ats = (nf.img);
BW_hpo = ((ats==rhp) | (ats==lhp)) ;
string_compcor = [darya_atlas_09c 'mni09c_hipp_mask.nii.gz'];
save_nii(make_nii(single(BW_hpo)),string_compcor)
unix(['fslcpgeom ' darya_atlas_09c 'mni_icbm152_CerebrA_tal_nlin_sym_09c.nii ' string_compcor])

%%create hippo masks
unix(['applywarp --ref=' dir_epi1 'func_mc_dc_tmean --in=' darya_atlas_09c 'mni09c_hipp_mask --out=' dir_epi1 'mni09c_hipp_2fmri --warp=' dir_epi1 'func2atlas_warp_inv --interp=nn'])
unix(['applywarp --ref=' dir_epi2 'func_mc_dc_tmean --in=' darya_atlas_09c 'mni09c_hipp_mask --out=' dir_epi2 'mni09c_hipp_2fmri --warp=' dir_epi2 'func2atlas_warp_inv --interp=nn'])
unix(['applywarp --ref=' dir_epi3 'func_mc_dc_tmeanc --in=' darya_atlas_09c 'mni09c_hipp_mask --out=' dir_epi3 'mni09c_hipp_2fmri --warp=' dir_epi3 'func2atlas_warp_inv --interp=nn'])

%% tsnr calculation for hippocampus
string_epi= 'func_mc_dc_hp_snr.nii.gz'; 
a0m=dir([dir_epi1,string_epi]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
epi1_dc = (nf.img);

string_epi= 'func_mc_dc_hp_physio_snr.nii.gz'; 
a0m=dir([dir_epi1,string_epi]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
epi1_dc_physio = (nf.img);

string_epi= 'func_mc_dc_hp_snr.nii.gz';
a0m=dir([dir_epi2,string_epi]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
epi2_dc = (nf.img);

string_epi= 'func_mc_dc_hp_physio_snr.nii.gz'; 
a0m=dir([dir_epi2,string_epi]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
epi2_dc_physio = (nf.img);

string_epi= 'func_mc_dc_hp_snr.nii.gz'; 
a0m=dir([dir_epi3,string_epi]);
nf=load_untouch_nii([dir_epi3,a0m(1).name]);
epi3_dc = (nf.img);

string_epi= 'func_mc_dc_hpc_physio_snr.nii.gz'; 
a0m=dir([dir_epi3,string_epi]);
nf=load_untouch_nii([dir_epi3,a0m(1).name]);
epi3_dc_physio = (nf.img);
%%rois
string_ats= 'mni09c_hipp_2fmri.nii.gz'; %% atlas in functional space for parcellation
a0m=dir([dir_epi1,string_ats]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
ats = (nf.img);
BW_hp_ap1 = (ats==1);

string_ats= 'mni09c_hipp_2fmri.nii.gz'; %% atlas in functional space for parcellation
a0m=dir([dir_epi2,string_ats]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
ats = (nf.img);
BW_hp_ap2 = (ats==1);

string_ats= 'mni09c_hipp_2fmri.nii.gz'; %% atlas in functional space for parcellation
a0m=dir([dir_epi3,string_ats]);
nf=load_untouch_nii([dir_epi3,a0m(1).name]);
ats = (nf.img);
BW_hp_ap3 = (ats==1);

epi1_dc_hp =mean(epi1_dc(BW_hp_ap1==1));
epi1_dc_physio_hp =mean(epi1_dc_physio(BW_hp_ap1==1));
epi2_dc_hp =mean(epi2_dc(BW_hp_ap2==1));
epi2_dc_physio_hp =mean(epi2_dc_physio(BW_hp_ap2==1));
epi3_dc_hp =mean(epi3_dc(BW_hp_ap3==1));
epi3_dc_physio_hp =mean(epi3_dc_physio(BW_hp_ap3==1));


hp_tsnr_rois = [epi1_dc_hp  epi1_dc_physio_hp;  epi2_dc_hp  epi2_dc_physio_hp; epi3_dc_hp  epi3_dc_physio_hp]

%% MNI_2009c Entorhinal mask labels

rec = 36;
lec = 87;
string_ats= 'mni_icbm152_CerebrA_tal_nlin_sym_09c.nii'; 
a0m=dir([darya_atlas_09c,string_ats]);
nf=load_untouch_nii([darya_atlas_09c,a0m(1).name]);
ats = (nf.img);
BW_ec = ((ats==rec) | (ats==lec)) ;
string_compcor = [darya_atlas_09c 'mni09c_ec_mask.nii.gz'];
save_nii(make_nii(single(BW_ec)),string_compcor)
unix(['fslcpgeom ' darya_atlas_09c 'mni_icbm152_CerebrA_tal_nlin_sym_09c.nii ' string_compcor])

%%create entorhinal masks
unix(['applywarp --ref=' dir_epi1 'func_mc_dc_tmean --in=' darya_atlas_09c 'mni09c_ec_mask --out=' dir_epi1 'mni09c_ec_2fmri --warp=' dir_epi1 'func2atlas_warp_inv --interp=nn'])
unix(['applywarp --ref=' dir_epi2 'func_mc_dc_tmean --in=' darya_atlas_09c 'mni09c_ec_mask --out=' dir_epi2 'mni09c_ec_2fmri --warp=' dir_epi2 'func2atlas_warp_inv --interp=nn'])
unix(['applywarp --ref=' dir_epi3 'func_mc_dc_tmeanc --in=' darya_atlas_09c 'mni09c_ec_mask --out=' dir_epi3 'mni09c_ec_2fmri --warp=' dir_epi3 'func2atlas_warp_inv --interp=nn'])


%%rois
string_ats= 'mni09c_ec_2fmri.nii.gz'; %% atlas in functional space for parcellation
a0m=dir([dir_epi1,string_ats]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
ats = (nf.img);
BW_ec_ap1 = (ats==1);

string_ats= 'mni09c_ec_2fmri.nii.gz'; %% atlas in functional space for parcellation
a0m=dir([dir_epi2,string_ats]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
ats = (nf.img);
BW_ec_ap2 = (ats==1);

string_ats= 'mni09c_ec_2fmri.nii.gz'; %% atlas in functional space for parcellation
a0m=dir([dir_epi3,string_ats]);
nf=load_untouch_nii([dir_epi3,a0m(1).name]);
ats = (nf.img);
BW_ec_ap3 = (ats==1);

epi1_dc_ec =mean(epi1_dc(BW_ec_ap1==1));
epi1_dc_physio_ec =mean(epi1_dc_physio(BW_ec_ap1==1));
epi2_dc_ec =mean(epi2_dc(BW_ec_ap2==1));
epi2_dc_physio_ec =mean(epi2_dc_physio(BW_ec_ap2==1));
epi3_dc_ec =mean(epi3_dc(BW_ec_ap3==1));
epi3_dc_physio_ec =mean(epi3_dc_physio(BW_ec_ap3==1));


hp_tsnr_rois = [epi1_dc_ec  epi1_dc_physio_ec;  epi2_dc_ec  epi2_dc_physio_ec; epi3_dc_ec  epi3_dc_physio_ec]

string_compcor = [dir_epi1 'func_hipp_mask.nii.gz'];
save_nii(make_nii(single(BW_hp_ap)),string_compcor)
unix(['fslcpgeom ' dir_epi1 'AP_mc_dc_tmean.nii.gz ' string_compcor])

string_compcor = [dir_epi2 'func_hipp_mask.nii.gz'];
save_nii(make_nii(single(BW_hp_pa)),string_compcor)
unix(['fslcpgeom ' dir_epi2 'PA_mc_dc_tmean.nii.gz ' string_compcor])

string_compcor = [dir_epi2 'func_bs_mask.nii.gz'];
save_nii(make_nii(single(BW_st_pa)),string_compcor)
unix(['fslcpgeom ' dir_epi2 'PA_mc_dc_tmean.nii.gz ' string_compcor])

%% functional analysis & registeration

%% task regressor in the fsl format
task_dir = '/export/home/vmalekian/LC_fMRI/misun/task/';
run1 = load([task_dir 'glmVer3_condi_run3.mat'])
dur1= run1.durations{1};
dur2 = run1.durations{2};
ons1= run1.onsets{1};
ons2 = run1.onsets{2};
zer1 = ones(length(ons1),1);
zer2 = ones(length(ons2),1);

cont1 = cat(2,ons1,dur1,zer1);
cont2 = cat(2,ons2,dur2,zer2);

filename = [task_dir 'run3_cont1.txt'];
fileID = fopen(filename,'w');
for i = 1:size(cont1, 1)
    fprintf(fileID, '%f ', cont1(i, :));
    fprintf(fileID, '\n');
end
fclose(fileID);

filename = [task_dir 'run3_cont2.txt'];
fileID = fopen(filename,'w');
for i = 1:size(cont2, 1)
    fprintf(fileID, '%f ', cont2(i, :));
    fprintf(fileID, '\n');
end
fclose(fileID);

%%run staticital analysis using fsl FEAT GUI
%%register functional maps to M2PRAGE space 

dir_feat1 = [ dir_epi1 'main.feat/'];
dir_feat2 = [ dir_epi2 'main.feat/'];
dir_feat3 = [ dir_epi3 'main.feat/'];
dir_feat1p = [ dir_epi1 'main_physio.feat/'];
dir_feat2p = [ dir_epi2 'main_physio.feat/'];
dir_feat3p = [ dir_epi3 'main_physio.feat/'];

unix(['flirt -in ' dir_feat1 'thresh_zstat1 -ref ' dir_struct 'T1_reo.anat/T1_biascorr_brain -out ' dir_feat1 'thresh_zstats1_highres -applyxfm -init ' dir_epi1 'func_mc_dc_tmean_highres.mat -interp spline'])
unix(['flirt -in ' dir_feat2 'thresh_zstat1 -ref ' dir_struct 'T1_reo.anat/T1_biascorr_brain -out ' dir_feat2 'thresh_zstats1_highres -applyxfm -init ' dir_epi2 'func_mc_dc_tmean_highres.mat -interp spline'])
unix(['flirt -in ' dir_feat3 'thresh_zstat1 -ref ' dir_struct 'T1_reo.anat/T1_biascorr_brain -out ' dir_feat3 'thresh_zstats1_highres -applyxfm -init ' dir_epi3 'func_mc_dc_tmean_highres.mat -interp spline'])

unix(['flirt -in ' dir_feat1p 'thresh_zstat1 -ref ' dir_struct 'T1_reo.anat/T1_biascorr_brain -out ' dir_feat1p 'thresh_zstats1_highres -applyxfm -init ' dir_epi1 'func_mc_dc_tmean_highres.mat -interp spline'])
unix(['flirt -in ' dir_feat2p 'thresh_zstat1 -ref ' dir_struct 'T1_reo.anat/T1_biascorr_brain -out ' dir_feat2p 'thresh_zstats1_highres -applyxfm -init ' dir_epi2 'func_mc_dc_tmean_highres.mat -interp spline'])
unix(['flirt -in ' dir_feat3p 'thresh_zstat1 -ref ' dir_struct 'T1_reo.anat/T1_biascorr_brain -out ' dir_feat3p 'thresh_zstats1_highres -applyxfm -init ' dir_epi3 'func_mc_dc_tmean_highres.mat -interp spline'])

unix(['flirt -in ' dir_feat1 'thresh_zstat2 -ref ' dir_struct 'T1_reo.anat/T1_biascorr_brain -out ' dir_feat1 'thresh_zstats2_highres -applyxfm -init ' dir_epi1 'func_mc_dc_tmean_highres.mat -interp spline'])
unix(['flirt -in ' dir_feat2 'thresh_zstat2 -ref ' dir_struct 'T1_reo.anat/T1_biascorr_brain -out ' dir_feat2 'thresh_zstats2_highres -applyxfm -init ' dir_epi2 'func_mc_dc_tmean_highres.mat -interp spline'])
unix(['flirt -in ' dir_feat3 'thresh_zstat2 -ref ' dir_struct 'T1_reo.anat/T1_biascorr_brain -out ' dir_feat3 'thresh_zstats2_highres -applyxfm -init ' dir_epi3 'func_mc_dc_tmean_highres.mat -interp spline'])

unix(['flirt -in ' dir_feat1p 'thresh_zstat2 -ref ' dir_struct 'T1_reo.anat/T1_biascorr_brain -out ' dir_feat1p 'thresh_zstats2_highres -applyxfm -init ' dir_epi1 'func_mc_dc_tmean_highres.mat -interp spline'])
unix(['flirt -in ' dir_feat2p 'thresh_zstat2 -ref ' dir_struct 'T1_reo.anat/T1_biascorr_brain -out ' dir_feat2p 'thresh_zstats2_highres -applyxfm -init ' dir_epi2 'func_mc_dc_tmean_highres.mat -interp spline'])
unix(['flirt -in ' dir_feat3p 'thresh_zstat2 -ref ' dir_struct 'T1_reo.anat/T1_biascorr_brain -out ' dir_feat3p 'thresh_zstats2_highres -applyxfm -init ' dir_epi3 'func_mc_dc_tmean_highres.mat -interp spline'])

unix(['applywarp --in=' dir_feat1p 'thresh_zstat2 --ref=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c_masked --out=' dir_feat1p 'thresh_zstat2_mni --warp=' dir_epi1 'func2atlas_warp --interp=spline'])
unix(['applywarp --in=' dir_feat2p 'thresh_zstat2 --ref=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c_masked --out=' dir_feat2p 'thresh_zstat2_mni --warp=' dir_epi2 'func2atlas_warp --interp=spline'])
unix(['applywarp --in=' dir_feat3p 'thresh_zstat2 --ref=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c_masked --out=' dir_feat3p 'thresh_zstat2_mni --warp=' dir_epi3 'func2atlas_warp --interp=spline'])
%%



