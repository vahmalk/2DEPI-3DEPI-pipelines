%% Pipeline for assessing 2D EPI data in Hippocampus,LC and Brainstem
% - AP vs PA phase encoding direction
% - tSNR, local TE and Jacobian
% - no cleaning vs cleaning via aCompCor and Spike
% - using Draya's Frank data 
% V Malekian, FIL Physics, 03/06/2024

clc
close all
clear all

addpath('/export/home/vmalekian/Malab_lib/NIfTI_20140122/')
addpath(genpath('/export/home/vmalekian/Yan/physio_correction/fmri_denoising-master'))

%% MP2RAGE preparations and segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_inv2 = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_largerFOV_INV2_0006/';
dir_struct = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_largerFOV_UNI_Images_0008/';
unix(['fslreorient2std ' dir_struct 's*1.nii ' dir_struct 'T1_reo.nii.gz'])
unix(['fslreorient2std ' dir_inv2 's*1.nii ' dir_inv2 'T1_reo.nii.gz'])

unix(['fsl_anat -i ' dir_struct 'T1_reo -o ' dir_struct 'T1_reo --noreorient --nocrop'])
unix(['fsl_anat -i ' dir_inv2 'T1_reo -o ' dir_inv2 'T1_reo --noreorient --nocrop'])

unix(['fslmaths ' dir_struct 'T1_reo.anat/T1_biascorr -mul ' dir_inv2 'T1_reo.anat/T1_biascorr_brain_mask ' dir_struct 'T1_reo.anat/T1_biascorr_brain'])


%% SPM structural tissue segmentation
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
unix(['fsl_anat -i ' dir_struct 'T1_reo_brain -o ' dir_struct 'T1_reo_brain --noreorient --nocrop'])


%% B0 field mapping preparation (phase and magnitude)

dir_mg = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/gre_1441B_B0_40mm_0002/';
dir_ph = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/gre_1441B_B0_40mmB0_Map_0003/';

unix(['fslreorient2std ' dir_mg 's*1.nii ' dir_mg 'magn1_reo.nii.gz'])
unix(['fslreorient2std ' dir_mg 's*2.nii ' dir_mg 'magn2_reo.nii.gz'])
unix(['fslmaths ' dir_mg 'magn1_reo -add ' dir_mg 'magn2_reo -div 2 ' dir_mg 'magn_reo'])
unix(['fslreorient2std ' dir_ph 's*.nii ' dir_ph 'phase_reo.nii.gz'])

unix(['bet ' dir_mg 'magn_reo.nii.gz ' dir_mg 'magn_reo_brain1.nii.gz -B -f 0.5 -g 0 -m'])
unix(['fslmaths ' dir_mg 'magn_reo_brain1.nii.gz -ero ' dir_mg 'magn_reo_brain.nii.gz'])
unix(['fsl_prepare_fieldmap SIEMENS ' dir_ph 'phase_reo.nii.gz ' dir_mg 'magn_reo_brain.nii.gz ' dir_ph 'fmap_rads.nii.gz 3.08'])
unix(['fslmaths ' dir_ph 'fmap_rads.nii.gz -div 6.28 ' dir_ph 'fmap_hz.nii.gz'])

unix(['fslmaths ' dir_epi1 'func_mc -Tmean ' dir_epi1 'func_mc_tmean'])
unix(['epi_reg --epi=' dir_epi1 'func_mc_tmean --t1=' dir_struct 'T1_reo.anat/T1_biascorr --t1brain=' dir_struct 'T1_reo.anat/T1_biascorr_brain --out=' dir_epi1 'func_mc_tmean_highres' ...
    ' --fmap=' dir_ph 'fmap_rads --fmapmag=' dir_mg 'magn_reo --fmapmagbrain=' dir_mg 'magn_reo_brain --echospacing=0.000416 --pedir=-y'])
unix(['fslmaths ' dir_epi2 'func_mc -Tmean ' dir_epi2 'func_mc_tmean'])
unix(['epi_reg --epi=' dir_epi2 'func_mc_tmean --t1=' dir_struct 'T1_reo.anat/T1_biascorr --t1brain=' dir_struct 'T1_reo.anat/T1_biascorr_brain --out=' dir_epi2 'func_mc_tmean_highres' ...
    ' --fmap=' dir_ph 'fmap_rads --fmapmag=' dir_mg 'magn_reo --fmapmagbrain=' dir_mg 'magn_reo_brain --echospacing=0.000416 --pedir=y'])

unix(['bet ' dir_epi1 'func_mc_tmean ' dir_epi1 'func_mc_tmean_brain -f 0.25 -g 0 -m'])
unix(['bet ' dir_epi2 'func_mc_tmean ' dir_epi2 'func_mc_tmean_brain -f 0.25 -g 0 -m'])

unix(['fslmaths ' dir_epi1 'func_mc_tmean_highres_fieldmaprads2epi.nii.gz -mas ' dir_epi1 'func_mc_tmean_brain_mask ' dir_epi1 'func_mc_tmean_highres_fieldmaprads2epi_masked'])
unix(['fslmaths ' dir_epi2 'func_mc_tmean_highres_fieldmaprads2epi.nii.gz -mas ' dir_epi2 'func_mc_tmean_brain_mask ' dir_epi2 'func_mc_tmean_highres_fieldmaprads2epi_masked'])

unix(['fslmaths ' dir_epi1 'func_mc_tmean_highres_fieldmaprads2epi_masked.nii.gz -div 6.28 ' dir_epi1 'func_mc_tmean_highres_fieldmaprads2epi_masked_hz.nii.gz'])
unix(['fslmaths ' dir_epi2 'func_mc_tmean_highres_fieldmaprads2epi_masked.nii.gz -div 6.28 ' dir_epi2 'func_mc_tmean_highres_fieldmaprads2epi_masked_hz.nii.gz'])


%% Start of the pipeline for 2D-EPI AP & PA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_epi1 = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/cmrr_mbep2d_1point25_hiedi_GRE_108_AP_0014/';
dir_epi2 = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/cmrr_mbep2d_1point25_hiedi_GRE_108_PA_0016/';

dir_distortion = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/distortion_rF_1/';

unix(['fslmerge -t ' dir_epi1 'func ' dir_epi1 'f*.nii'])
unix(['fslmerge -tr ' dir_epi1 'func ' dir_epi1 'func 2.1'])

unix(['fslmerge -t ' dir_epi2 'func ' dir_epi2 'f*.nii'])
unix(['fslmerge -tr ' dir_epi2 'func ' dir_epi2 'func 2.1'])

unix(['mcflirt -in ' dir_epi1 'func -o ' dir_epi1 'func_mc -refvol 219 -spline_final -plots -report -stats -mats'])
unix(['mcflirt -in ' dir_epi2 'func -o ' dir_epi2 'func_mc -refvol 0 -spline_final -plots -report -stats -mats'])
unix(['fslmaths ' dir_epi1 'func_mc -thr 0 ' dir_epi1 'func_mc'])
unix(['fslmaths ' dir_epi2 'func_mc -thr 0 ' dir_epi2 'func_mc'])

unix(['fslroi ' dir_epi1 'func_mc ' dir_epi1 'AP_lst 0 -1 0 -1 0 -1 219 1'])
unix(['fslroi ' dir_epi2 'func_mc ' dir_epi2 'PA_fst 0 -1 0 -1 0 -1 0 1'])

%%for AP
unix(['fslmerge -t ' dir_distortion 'AP_PA ' dir_epi1 'AP_lst ' dir_epi2 'PA_fst'])
dir_distortion1 = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/distortion_rF_1/';
unix(['topup --imain=' dir_distortion 'AP_PA --datain=' dir_distortion1 'param.txt --config=b02b0.cnf --out=' dir_distortion1 'topup_AP_PA --fout=' dir_distortion1 'topup_field ' ...
    '--iout=' dir_distortion1 'topup_unwarp --jacout=' dir_distortion1 'jac --rbmout=' dir_distortion1 'xfm --dfout=' dir_distortion1 'topup_warpfield'])
unix(['applytopup --imain=' dir_epi1 'func_mc --topup=' dir_distortion1 'topup_AP_PA --datain=' dir_distortion1 'param.txt --inindex=1 --out=' dir_epi1 'AP_mc_dc --method=jac']) %% applytopup with Jacobian correction
unix(['fslmaths ' dir_epi1 'AP_mc_dc -thr 0 ' dir_epi1 'func_mc_dc'])
%%for PA
unix(['fslmerge -t ' dir_distortion 'PA_AP ' dir_epi2 'PA_fst ' dir_epi1 'AP_lst'])
dir_distortion2 = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/distortion_rF_2/';
unix(['topup --imain=' dir_distortion 'PA_AP --datain=' dir_distortion2 'param_PA.txt --config=b02b0.cnf --out=' dir_distortion2 'topup_PA_AP --fout=' dir_distortion2 'topup_field_PA ' ...
    '--iout=' dir_distortion2 'topup_unwarp_PA --jacout=' dir_distortion2 'jac_PA --rbmout=' dir_distortion2 'xfm --dfout=' dir_distortion2 'topup_warpfield_PA'])
unix(['applytopup --imain=' dir_epi2 'func_mc --topup=' dir_distortion2 'topup_PA_AP --datain=' dir_distortion2 'param_PA.txt --inindex=1 --out=' dir_epi2 'PA_mc_dc --method=jac']) %% applytopup with Jacobian correction
unix(['fslmaths ' dir_epi2 'PA_mc_dc -thr 0 ' dir_epi2 'func_mc_dc'])


%% Perform HighPass, BET, Co-registeration, Normalisation
unix(['fslmaths ' dir_epi1 'AP_mc_dc -Tmean ' dir_epi1 'AP_mc_dc_tmean'])
unix(['fslmaths ' dir_epi2 'PA_mc_dc -Tmean ' dir_epi2 'PA_mc_dc_tmean'])

unix (['fslmaths ' dir_epi1 'AP_mc_dc -bptf 23.809 -1 -add ' dir_epi1 'AP_mc_dc_tmean ' dir_epi1 'AP_mc_dc_hp']) % 1/(2*0.01*2.1)
unix (['fslmaths ' dir_epi2 'PA_mc_dc -bptf 23.809 -1 -add ' dir_epi2 'PA_mc_dc_tmean ' dir_epi2 'PA_mc_dc_hp']) % 1/(2*0.01*2.1)

unix(['fslmaths ' dir_epi1 'AP_mc_dc_hp -thr 0 ' dir_epi1 'AP_mc_dc_hp'])
unix(['fslmaths ' dir_epi2 'PA_mc_dc_hp -thr 0 ' dir_epi2 'PA_mc_dc_hp'])


%% Structral registeration BBR

%%without B0 field map
unix(['epi_reg --epi=' dir_epi1 'AP_mc_dc_tmean --t1=' dir_struct 'T1_reo.anat/T1_biascorr --t1brain=' dir_struct 'T1_reo.anat/T1_biascorr_brain --out=' dir_epi1 'AP_mc_dc_tmean_highres --pedir=y-'])
unix(['epi_reg --epi=' dir_epi2 'PA_mc_dc_tmean --t1=' dir_struct 'T1_reo.anat/T1_biascorr --t1brain=' dir_struct 'T1_reo.anat/T1_biascorr_brain --out=' dir_epi2 'PA_mc_dc_tmean_highres --pedir=y'])

%%with B0 field map
unix(['epi_reg --epi=' dir_epi1 'AP_cut_mc_tmean --t1=' dir_struct 'T1.anat/T1_biascorr --t1brain=' dir_struct 'T1.anat/T1_biascorr_brain --out=' dir_epi1 'AP_cut_mc_tmean_highres' ...
    ' --fmap=' dir_ph 'fmap_rads --fmapmag=' dir_mg 'magn_reo --fmapmagbrain=' dir_mg 'magn_reo_brain --echospacing=0.000365 --pedir=-y'])
unix(['epi_reg --epi=' dir_epi2 'PA_cut_mc_tmean --t1=' dir_struct 'T1.anat/T1_biascorr --t1brain=' dir_struct 'T1.anat/T1_biascorr_brain --out=' dir_epi2 'PA_cut_mc_final_tmean_highres' ...
    ' --fmap=' dir_ph 'fmap_rads --fmapmag=' dir_mg 'magn_reo --fmapmagbrain=' dir_mg 'magn_reo_brain --echospacing=0.000365 --pedir=y'])


%% aCompCor for extracting physio regressers from CSF and WM
unix(['fslmaths ' dir_struct 'T1_reo.anat/c2T1_biascorr -thr 0.99 -bin ' dir_struct 'c2T1_reo_mask']) %% WM
unix(['fslmaths ' dir_struct 'T1_reo.anat/c3T1_biascorr -thr 0.99 -bin ' dir_struct 'c3T1_reo_mask']) %% CSF

unix(['bet ' dir_epi1 'AP_mc_dc_tmean ' dir_epi1 'AP_mc_dc_tmean_brain -f 0.5 -g 0 -m'])
unix(['bet ' dir_epi2 'PA_mc_dc_tmean ' dir_epi2 'PA_mc_dc_tmean_brain -f 0.5 -g 0 -m'])

unix(['convert_xfm -inverse -omat ' dir_epi1 'AP_mc_dc_tmean_highres_inv.mat ' dir_epi1 'AP_mc_dc_tmean_highres.mat'])
unix(['convert_xfm -inverse -omat ' dir_epi2 'PA_mc_dc_tmean_highres_inv.mat ' dir_epi2 'PA_mc_dc_tmean_highres.mat'])

unix(['applywarp --ref=' dir_epi1 'AP_mc_dc_tmean --in=' dir_struct 'c2T1_reo_mask' ...
    ' --out=' dir_epi1 'c2T1_reo_2fMREPI --premat=' dir_epi1 'AP_mc_dc_tmean_highres_inv.mat --interp=nn'])
unix(['applywarp --ref=' dir_epi1 'AP_mc_dc_tmean --in=' dir_struct 'c3T1_reo_mask' ...
    ' --out=' dir_epi1 'c3T1_reo_2fMREPI --premat=' dir_epi1 'AP_mc_dc_tmean_highres_inv.mat --interp=nn'])

unix(['applywarp --ref=' dir_epi2 'PA_mc_dc_tmean --in=' dir_struct 'c2T1_reo_mask' ...
    ' --out=' dir_epi2 'c2T1_reo_2fMREPI --premat=' dir_epi2 'PA_mc_dc_tmean_highres_inv.mat --interp=nn'])
unix(['applywarp --ref=' dir_epi2 'PA_mc_dc_tmean --in=' dir_struct 'c3T1_reo_mask' ...
    ' --out=' dir_epi2 'c3T1_reo_2fMREPI --premat=' dir_epi2 'PA_mc_dc_tmean_highres_inv.mat --interp=nn'])

unix(['fslmaths ' dir_epi1 'c2T1_reo_2fMREPI -mul ' dir_epi1 'AP_mc_dc_tmean_brain_mask ' dir_epi1 'c2T1_reo_2fMREPI_BET'])
unix(['fslmaths ' dir_epi1 'c3T1_reo_2fMREPI -mul ' dir_epi1 'AP_mc_dc_tmean_brain_mask ' dir_epi1 'c3T1_reo_2fMREPI_BET'])

unix(['fslmaths ' dir_epi2 'c2T1_reo_2fMREPI -mul ' dir_epi2 'PA_mc_dc_tmean_brain_mask ' dir_epi2 'c2T1_reo_2fMREPI_BET'])
unix(['fslmaths ' dir_epi2 'c3T1_reo_2fMREPI -mul ' dir_epi2 'PA_mc_dc_tmean_brain_mask ' dir_epi2 'c3T1_reo_2fMREPI_BET'])

string = 'c2T1_reo_2fMREPI_BET.nii.gz';
a1m=dir([dir_epi1,string]);
nf=load_untouch_nii([dir_epi1,a1m(1).name]);
roi_W= (nf.img);
string = 'c3T1_reo_2fMREPI_BET.nii.gz';
a1m=dir([dir_epi1,string]);
nf=load_untouch_nii([dir_epi1,a1m(1).name]);
roi_C= (nf.img);

string = 'AP_mc_dc_hp.nii.gz';
a1m=dir([dir_epi1,string]);
nf=load_untouch_nii([dir_epi1,a1m(1).name]);
data= (nf.img);
rois1 = {roi_C,roi_W};
X= fmri_compcor(data,rois1,[7 7]);% Alternitvely can use X= fmri_compcor(data,rois,[0.25 0.25])
[res,model_info] = fmri_cleaning(data,-1,X);
string_compcor = [dir_epi1 'AP_mc_dc_hp_physio.nii.gz'];
save_nii(make_nii((res)),string_compcor)
unix(['fslcpgeom ' dir_epi1 'AP_mc_dc_hp.nii.gz ' string_compcor])

figure,
for k=1:14
 subplot(7,2,k),plot(X(:,k))
end

string = 'c2T1_reo_2fMREPI_BET.nii.gz';
a1m=dir([dir_epi2,string]);
nf=load_untouch_nii([dir_epi2,a1m(1).name]);
roi_W= (nf.img);
string = 'c3T1_reo_2fMREPI_BET.nii.gz';
a1m=dir([dir_epi2,string]);
nf=load_untouch_nii([dir_epi2,a1m(1).name]);
roi_C= (nf.img);
string = 'PA_mc_dc_hp.nii.gz';
a1m=dir([dir_epi2,string]);
nf=load_untouch_nii([dir_epi2,a1m(1).name]);
data= (nf.img);
rois2 = {roi_C,roi_W};
X= fmri_compcor(data,rois2,[7 7]);% Alternitvely can use X= fmri_compcor(data,rois,[0.25 0.25])
[res,model_info] = fmri_cleaning(data,-1,X);
string_compcor = [dir_epi2 'PA_mc_dc_hp_physio.nii.gz'];
save_nii(make_nii((res)),string_compcor)
unix(['fslcpgeom ' dir_epi2 'PA_mc_dc_hp.nii.gz ' string_compcor])

figure
for k=1:14
 subplot(7,2,k),plot(X(:,k))
end

%% Load physio regressers from physio matlab code to clean fMRI data (spike)
load('/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/physio/restricted_AP_run1_R_session1.mat')   % EPI AP
Rap = R;
string = 'AP_mc_dc_hp.nii.gz';
a1m=dir([dir_epi1,string]);
nf=load_untouch_nii([dir_epi1,a1m(1).name]);
data= (nf.img);
[res,model_info] = fmri_cleaning(data,-1,Rap);
string_compcor = [dir_epi1 'AP_mc_dc_hp_spike.nii.gz'];
save_nii(make_nii((res)),string_compcor)
unix(['fslcpgeom ' dir_epi1 'AP_mc_dc_hp.nii.gz ' string_compcor])

figure
for k=1:14
 subplot(7,2,k),plot(Rap(:,k))
end


load('/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/physio/restricted_PA_run2_R_session1.mat')   % EPI PA
Rpa = R;
string = 'PA_mc_dc_hp.nii.gz';
a1m=dir([dir_epi2,string]);
nf=load_untouch_nii([dir_epi2,a1m(1).name]);
data= (nf.img);
[res,model_info] = fmri_cleaning(data,-1,Rpa);
string_compcor = [dir_epi2 'PA_mc_dc_hp_spike.nii.gz'];
save_nii(make_nii((res)),string_compcor)
unix(['fslcpgeom ' dir_epi2 'PA_mc_dc_hp.nii.gz ' string_compcor])

figure
for k=1:14
 subplot(7,2,k),plot(Rpa(:,k))
end

%% tsnr before cleaning
unix(['fslmaths ' dir_epi1 'AP_mc_dc_hp -Tmean ' dir_epi1 'AP_mc_dc_hp_tmean'])
unix(['fslmaths ' dir_epi1 'AP_mc_dc_hp -Tstd ' dir_epi1 'AP_mc_dc_hp_std'])
unix(['fslmaths ' dir_epi1 'AP_mc_dc_hp_tmean -div ' dir_epi1 'AP_mc_dc_hp_std ' dir_epi1 'AP_mc_dc_hp_snr'])

unix(['fslmaths ' dir_epi2 'PA_mc_dc_hp -Tmean ' dir_epi2 'PA_mc_dc_hp_tmean'])
unix(['fslmaths ' dir_epi2 'PA_mc_dc_hp -Tstd ' dir_epi2 'PA_mc_dc_hp_std'])
unix(['fslmaths ' dir_epi2 'PA_mc_dc_hp_tmean -div ' dir_epi2 'PA_mc_dc_hp_std ' dir_epi2 'PA_mc_dc_hp_snr'])

%% tsnr after cleaning using acompcor

unix(['fslmaths ' dir_epi1 'AP_mc_dc_hp_physio -Tmean ' dir_epi1 'AP_mc_dc_hp_physio_tmean'])
unix(['fslmaths ' dir_epi1 'AP_mc_dc_hp_physio -Tstd ' dir_epi1 'AP_mc_dc_hp_physio_std'])
unix(['fslmaths ' dir_epi1 'AP_mc_dc_hp_physio_tmean -div ' dir_epi1 'AP_mc_dc_hp_physio_std ' dir_epi1 'AP_mc_dc_hp_physio_snr'])

unix(['fslmaths ' dir_epi2 'PA_mc_dc_hp_physio -Tmean ' dir_epi2 'PA_mc_dc_hp_physio_tmean'])
unix(['fslmaths ' dir_epi2 'PA_mc_dc_hp_physio -Tstd ' dir_epi2 'PA_mc_dc_hp_physio_std'])
unix(['fslmaths ' dir_epi2 'PA_mc_dc_hp_physio_tmean -div ' dir_epi2 'PA_mc_dc_hp_physio_std ' dir_epi2 'PA_mc_dc_hp_physio_snr'])

%% tsnr after cleaning using spike

unix(['fslmaths ' dir_epi1 'AP_mc_dc_hp_spike -Tmean ' dir_epi1 'AP_mc_dc_hp_spike_tmean'])
unix(['fslmaths ' dir_epi1 'AP_mc_dc_hp_spike -Tstd ' dir_epi1 'AP_mc_dc_hp_spike_std'])
unix(['fslmaths ' dir_epi1 'AP_mc_dc_hp_spike_tmean -div ' dir_epi1 'AP_mc_dc_hp_spike_std ' dir_epi1 'AP_mc_dc_hp_spike_snr'])

unix(['fslmaths ' dir_epi2 'PA_mc_dc_hp_spike -Tmean ' dir_epi2 'PA_mc_dc_hp_spike_tmean'])
unix(['fslmaths ' dir_epi2 'PA_mc_dc_hp_spike -Tstd ' dir_epi2 'PA_mc_dc_hp_spike_std'])
unix(['fslmaths ' dir_epi2 'PA_mc_dc_hp_spike_tmean -div ' dir_epi2 'PA_mc_dc_hp_spike_std ' dir_epi2 'PA_mc_dc_hp_spike_snr'])


%% ROI definition for LC and hippocampus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ICBM atlas normalisation to extract LC 
darya_atlas = '/export/home/vmalekian/LC_fMRI/20231207.M700809_FIL/darya_atlas/';
unix(['flirt -in ' dir_struct 'T1_reo_brain -ref ' darya_atlas 'mni_icbm152_t1_tal_nlin_asym_09b_hires_FSL_bbox_struc_brain_CSFin -out ' dir_struct 'highres2icbm -omat ' dir_struct 'highres2icbm.mat -cost corratio -dof 12 -interp spline'])
unix(['fnirt --in=' dir_struct 'T1_reo_brain --ref=' darya_atlas 'mni_icbm152_t1_tal_nlin_asym_09b_hires_FSL_bbox_struc_brain_CSFin --aff=' dir_struct 'highres2icbm.mat' ...
    ' --cout=' dir_struct 'highres2icbm_warp --iout=' dir_struct 'highres2icbm_out --fout=' dir_struct 'highres2icbm_field --interp=spline']) 

%%AP
unix(['convertwarp --ref=' darya_atlas 'mni_icbm152_t1_tal_nlin_asym_09b_hires_FSL_bbox_struc_brain_CSFin --premat=' dir_epi1 'AP_mc_dc_tmean_highres.mat' ...
    ' --warp1=' dir_struct 'highres2icbm_warp --out=' dir_epi1 'example_func2icbm_warp'])
unix(['invwarp -w ' dir_epi1 'example_func2icbm_warp -o ' dir_epi1 'example_func2icbm_warp_inv -r ' dir_epi1 'AP_mc_dc_tmean'])
unix(['fslmaths ' darya_atlas 'LCTMP_n53_5SD_prob0 -bin ' darya_atlas 'LCTMP_n53_5SD_prob0_bin'])
unix(['applywarp --ref=' dir_epi1 'AP_mc_dc_tmean --in=' darya_atlas 'LCTMP_n53_5SD_prob0_bin' ...
    ' --out=' dir_epi1 'LCTMP_n53_5SD_prob0_2fmri --warp=' dir_epi1 'example_func2icbm_warp_inv --interp=nn'])
%%PA
unix(['convertwarp --ref=' darya_atlas 'mni_icbm152_t1_tal_nlin_asym_09b_hires_FSL_bbox_struc_brain_CSFin --premat=' dir_epi2 'PA_mc_dc_tmean_highres.mat' ...
    ' --warp1=' dir_struct 'highres2icbm_warp --out=' dir_epi2 'example_func2icbm_warp'])
unix(['invwarp -w ' dir_epi2 'example_func2icbm_warp -o ' dir_epi2 'example_func2icbm_warp_inv -r ' dir_epi2 'PA_mc_dc_tmean'])
unix(['applywarp --ref=' dir_epi2 'PA_mc_dc_tmean --in=' darya_atlas 'LCTMP_n53_5SD_prob0_bin' ...
    ' --out=' dir_epi2 'LCTMP_n53_5SD_prob0_2fmri --warp=' dir_epi2 'example_func2icbm_warp_inv --interp=nn'])

%% MNI_2009C atlas normalisation to extract brain stem & hippo

darya_atlas_09c = '/export/home/vmalekian/LC_fMRI/20231207.M700809_FIL/darya_atlas/mni_icbm152_nlin_asym_09c/';
unix(['flirt -in ' dir_struct 'T1_reo_brain -ref ' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c -out ' dir_struct 'highres2mni9c -omat ' dir_struct 'highres2mni9c.mat -cost corratio -dof 12 -interp spline'])
unix(['fnirt --in=' dir_struct 'T1_reo_brain --ref=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c --aff=' dir_struct 'highres2mni9c.mat' ...
    ' --cout=' dir_struct 'highres2mni9c_warp --iout=' dir_struct 'highres2mni9c_out --fout=' dir_struct 'highres2mni9c_field --interp=spline']) 
%%AP
unix(['convertwarp --ref=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c --premat=' dir_epi1 'AP_mc_dc_tmean_highres.mat' ...
    ' --warp1=' dir_struct 'highres2mni9c_warp --out=' dir_epi1 'example_func2mni09c_warp'])
unix(['invwarp -w ' dir_epi1 'example_func2mni09c_warp -o ' dir_epi1 'example_func2mni09c_warp_inv -r ' dir_epi1 'AP_mc_dc_tmean'])
unix(['fslmaths ' darya_atlas_09c 'SN_VTA_MNI2009c -bin ' darya_atlas_09c 'SN_VTA_MNI2009c_bin'])
unix(['applywarp --ref=' dir_epi1 'AP_mc_dc_tmean --in=' darya_atlas_09c 'SN_VTA_MNI2009c_bin' ...
    ' --out=' dir_epi1 'SN_VTA_MNI2009c_2fmri --warp=' dir_epi1 'example_func2mni09c_warp_inv --interp=nn'])
%%PA
unix(['convertwarp --ref=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c --premat=' dir_epi2 'PA_mc_dc_tmean_highres.mat' ...
    ' --warp1=' dir_struct 'highres2mni9c_warp --out=' dir_epi2 'example_func2mni09c_warp'])
unix(['invwarp -w ' dir_epi2 'example_func2mni09c_warp -o ' dir_epi2 'example_func2mni09c_warp_inv -r ' dir_epi2 'PA_mc_dc_tmean'])
unix(['applywarp --ref=' dir_epi2 'PA_mc_dc_tmean --in=' darya_atlas_09c 'SN_VTA_MNI2009c_bin' ...
    ' --out=' dir_epi2 'SN_VTA_MNI2009c_2fmri --warp=' dir_epi2 'example_func2mni09c_warp_inv --interp=nn'])

unix(['applywarp --ref=' darya_atlas_09c 'mni_icbm152_t1_tal_nlin_asym_09c --in=' dir_epi2 'PA_mc_dc_hp' ...
    ' --out=' dir_epi2 'PA_mc_dc_hp2mni09c --warp=' dir_epi2 'example_func2mni09c_warp --interp=spline'])
unix(['fslmaths ' dir_epi2 'PA_mc_dc_hp2mni09c -Tmean ' dir_epi2 'PA_mc_dc_hp2mni09c_tmean'])
unix(['fslmaths ' dir_epi2 'PA_mc_dc_hp2mni09c -Tstd ' dir_epi2 'PA_mc_dc_hp2mni09c_std'])
unix(['fslmaths ' dir_epi2 'PA_mc_dc_hp2mni09c_tmean -div ' dir_epi2 'PA_mc_dc_hp2mni09c_std ' dir_epi2 'PA_mc_dc_hp2mni09c_snr'])


%% MNI_2009C Hippo labels

rhp = 48;
lhp = 99;
string_ats= 'mni_icbm152_CerebrA_tal_nlin_sym_09c.nii'; %% atlas in MP2RAGE space for parcellation
a0m=dir([darya_atlas_09c,string_ats]);
nf=load_untouch_nii([darya_atlas_09c,a0m(1).name]);
ats = (nf.img);
BW_hpo = ((ats==rhp) | (ats==lhp)) ;
string_compcor = [darya_atlas_09c 'mni09c_hipp_mask.nii.gz'];
save_nii(make_nii(single(BW_hpo)),string_compcor)
unix(['fslcpgeom ' darya_atlas_09c 'mni_icbm152_CerebrA_tal_nlin_sym_09c.nii ' string_compcor])

%%AP
unix(['applywarp --ref=' dir_epi1 'AP_mc_dc_tmean --in=' darya_atlas_09c 'mni09c_hipp_mask' ...
    ' --out=' dir_epi1 'mni09c_hipp_2fmri --warp=' dir_epi1 'example_func2mni09c_warp_inv --interp=nn'])
%%PA
unix(['applywarp --ref=' dir_epi2 'PA_mc_dc_tmean --in=' darya_atlas_09c 'mni09c_hipp_mask' ...
    ' --out=' dir_epi2 'mni09c_hipp_2fmri --warp=' dir_epi2 'example_func2mni09c_warp_inv --interp=nn'])

%% AP direction ROI tsnr calculation

string_epi= 'AP_mc_dc_hp_tmean.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi1,string_epi]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
epi_dc = (nf.img);

string_epi= 'AP_mc_dc_hp_physio_tmean.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi1,string_epi]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
epi_dc_physio = (nf.img);

string_epi= 'AP_mc_dc_hp_spike_tmean.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi1,string_epi]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
epi_dc_spike = (nf.img);

%%rois
string_ats= 'mni09c_hipp_2fmri.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi1,string_ats]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
ats = (nf.img);
BW_hp_ap = (ats==1);

string_ats= 'SN_VTA_MNI2009c_2fmri.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi1,string_ats]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
ats = (nf.img);
BW_bs_ap = (ats==1);

string_lc = 'LCTMP_n53_5SD_prob0_2fmri.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi1,string_lc]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
ats = (nf.img);
BW_lc_ap = (ats==1);

epi_dc_stem =mean(epi_dc(BW_bs_ap==1));
epi_dc_physio_stem =mean(epi_dc_physio(BW_bs_ap==1));
epi_dc_spike_stem =mean(epi_dc_spike(BW_bs_ap==1));

epi_dc_hp =mean(epi_dc(BW_hp_ap==1));
epi_dc_physio_hp =mean(epi_dc_physio(BW_hp_ap==1));
epi_dc_spike_hp =mean(epi_dc_spike(BW_hp_ap==1));

epi_dc_lc =mean(epi_dc(BW_lc_ap==1));
epi_dc_physio_lc =mean(epi_dc_physio(BW_lc_ap==1));
epi_dc_spike_lc =mean(epi_dc_spike(BW_lc_ap==1));

AP_tsnr_rois = [epi_dc_stem  epi_dc_hp  epi_dc_lc; epi_dc_physio_stem epi_dc_physio_hp  epi_dc_physio_lc; epi_dc_spike_stem epi_dc_spike_hp  epi_dc_spike_lc]

%% PA direction ROI tsnr calculation
string_epi= 'PA_mc_dc_hp_tmean.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi2,string_epi]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
epi_dc = (nf.img);

string_epi= 'PA_mc_dc_hp_physio_tmean.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi2,string_epi]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
epi_dc_physio = (nf.img);

string_epi= 'PA_mc_dc_hp_spike_tmean.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi2,string_epi]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
epi_dc_spike = (nf.img);

%%rois
string_ats= 'mni09c_hipp_2fmri.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi2,string_ats]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
ats = (nf.img);
BW_hp_pa = (ats==1);

string_ats= 'SN_VTA_MNI2009c_2fmri.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi2,string_ats]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
ats = (nf.img);
BW_bs_pa = (ats==1);

string_lc = 'LCTMP_n53_5SD_prob0_2fmri.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi2,string_lc]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
ats = (nf.img);
BW_lc_pa = (ats==1);

epi_dc_stem =mean(epi_dc(BW_bs_pa==1));
epi_dc_physio_stem =mean(epi_dc_physio(BW_bs_pa==1));
epi_dc_spike_stem =mean(epi_dc_spike(BW_bs_pa==1));

epi_dc_hp =mean(epi_dc(BW_hp_pa==1));
epi_dc_physio_hp =mean(epi_dc_physio(BW_hp_pa==1));
epi_dc_spike_hp =mean(epi_dc_spike(BW_hp_pa==1));

epi_dc_lc =mean(epi_dc(BW_lc_pa==1));
epi_dc_physio_lc =mean(epi_dc_physio(BW_lc_pa==1));
epi_dc_spike_lc =mean(epi_dc_spike(BW_lc_pa==1));

PA_tsnr_rois = [epi_dc_stem  epi_dc_hp  epi_dc_lc; epi_dc_physio_stem epi_dc_physio_hp  epi_dc_physio_lc; epi_dc_spike_stem epi_dc_spike_hp  epi_dc_spike_lc]

%%save masks as nifti files
string_compcor = [dir_epi1 'func_hipp_mask.nii.gz'];
save_nii(make_nii(single(BW_hp_ap)),string_compcor)
unix(['fslcpgeom ' dir_epi1 'AP_mc_dc_tmean.nii.gz ' string_compcor])

string_compcor = [dir_epi2 'func_hipp_mask.nii.gz'];
save_nii(make_nii(single(BW_hp_pa)),string_compcor)
unix(['fslcpgeom ' dir_epi2 'PA_mc_dc_tmean.nii.gz ' string_compcor])

string_compcor = [dir_epi2 'func_bs_mask.nii.gz'];
save_nii(make_nii(single(BW_st_pa)),string_compcor)
unix(['fslcpgeom ' dir_epi2 'PA_mc_dc_tmean.nii.gz ' string_compcor])

%% Compute the Local-TE and Jacobian measures for each ROI using barbara's code

string_epi= 'TE_local.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi1,string_epi]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
TE_ap = (nf.img);

string_epi= 'Jacobian.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi1,string_epi]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
jac_ap = (nf.img);

TE_ap_lc =mean(TE_ap(BW_lc_ap==1));
jac_ap_lc =mean(jac_ap(BW_lc_ap==1));
TE_ap_bs =mean(TE_ap(BW_bs_ap==1));
jac_ap_bs =mean(jac_ap(BW_bs_ap==1));
TE_ap_hipo =mean(TE_ap(BW_hp_ap==1));
jac_ap_hipo =mean(jac_ap(BW_hp_ap==1));

AP_TE_jac = [TE_ap_lc jac_ap_lc; TE_ap_bs jac_ap_bs; TE_ap_hipo jac_ap_hipo]

string_epi= 'TE_local.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi2,string_epi]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
TE_pa = (nf.img);

string_epi= 'Jacobian.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi2,string_epi]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
jac_pa = (nf.img);

TE_pa_lc =mean(TE_pa(BW_lc_pa==1));
jac_pa_lc =mean(jac_pa(BW_lc_pa==1));
TE_pa_bs =mean(TE_pa(BW_bs_pa==1));
jac_pa_bs =mean(jac_pa(BW_bs_pa==1));
TE_pa_hipo =mean(TE_pa(BW_hp_pa==1));
jac_pa_hipo =mean(jac_pa(BW_hp_pa==1));

PA_TE_jac = [TE_pa_lc jac_pa_lc; TE_pa_bs jac_pa_bs; TE_pa_hipo jac_pa_hipo]


%% use fsl_anat segmentation to extract sub-cortical ROI using non-linear transformations

unix(['fsl_anat -i ' dir_struct '*.nii -o ' dir_struct 'T1'])
unix(['fsl_anat -i ' dir_inv2 '*.nii -o ' dir_inv2 'T1'])

fsl_anat_dir  = '/export/home/vmalekian/LC_fMRI/20231207.M700809_FIL/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_largerFOV_UNI_Images_0008/T1.anat/first_results/';
unix(['applywarp --ref=' dir_epi1 'AP_cut_mc_dc_tmean --in=' fsl_anat_dir 'T1_first_all_fast_firstseg' ...
    ' --out=' dir_epi1 'T1_first_all_fast_firstseg_2fMREPI --premat=' dir_epi1 'AP_cut_mc_dc_tmean_highres_inv.mat --interp=nn'])

unix(['applywarp --ref=' dir_epi2 'PA_cut_mc_final_dc_tmean --in=' fsl_anat_dir 'T1_first_all_fast_firstseg' ...
    ' --out=' dir_epi2 'T1_first_all_fast_firstseg_2fMREPI --premat=' dir_epi2 'PA_cut_mc_final_dc_tmean_highres_inv.mat --interp=nn'])

rhp = 53;lhp = 17;stem = 16;

string_ats= 'T1_first_all_fast_firstseg_2fMREPI.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi1,string_ats]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
ats = (nf.img);
BW_st = (ats==stem);
BW_rhp = (ats==rhp);
BW_lhp = (ats==lhp);

string_epi= 'AP_cut_mc_dc_hp_snr.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi1,string_epi]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
epi_dc = (nf.img);

string_epi= 'AP_cut_mc_dc_hp_physio_snr.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi1,string_epi]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
epi_dc_physio = (nf.img);

string_epi= 'AP_cut_mc_dc_hp_spike_snr.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi1,string_epi]);
nf=load_untouch_nii([dir_epi1,a0m(1).name]);
epi_dc_spike = (nf.img);


epi_dc_stem =mean(epi_dc(BW_st==1));
epi_dc_physio_stem =mean(epi_dc_physio(BW_st==1));
epi_dc_spike_stem =mean(epi_dc_spike(BW_st==1));

epi_dc_rhp =mean(epi_dc(BW_rhp==1));
epi_dc_physio_rhp =mean(epi_dc_physio(BW_rhp==1));
epi_dc_spike_rhp =mean(epi_dc_spike(BW_rhp==1));

epi_dc_lhp =mean(epi_dc(BW_lhp==1));
epi_dc_physio_lhp =mean(epi_dc_physio(BW_lhp==1));
epi_dc_spike_lhp =mean(epi_dc_spike(BW_lhp==1));


brain_subcortical_AP = [epi_dc_stem  epi_dc_rhp  epi_dc_lhp; epi_dc_physio_stem epi_dc_physio_rhp epi_dc_physio_lhp;epi_dc_spike_stem epi_dc_spike_rhp epi_dc_spike_lhp]

string_epi= 'PA_cut_mc_dc_hp_snr.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi2,string_epi]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
epi_dc = (nf.img);

string_epi= 'PA_cut_mc_dc_hp_physio_snr.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi2,string_epi]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
epi_dc_physio = (nf.img);

string_epi= 'PA_cut_mc_dc_hp_spike_snr.nii.gz'; %% atlas in MP2RAGE space for parcellation
a0m=dir([dir_epi2,string_epi]);
nf=load_untouch_nii([dir_epi2,a0m(1).name]);
epi_dc_spike = (nf.img);

epi_dc_stem =mean(epi_dc(BW_st==1));
epi_dc_physio_stem =mean(epi_dc_physio(BW_st==1));
epi_dc_spike_stem =mean(epi_dc_spike(BW_st==1));

epi_dc_rhp =mean(epi_dc(BW_rhp==1));
epi_dc_physio_rhp =mean(epi_dc_physio(BW_rhp==1));
epi_dc_spike_rhp =mean(epi_dc_spike(BW_rhp==1));

epi_dc_lhp =mean(epi_dc(BW_lhp==1));
epi_dc_physio_lhp =mean(epi_dc_physio(BW_lhp==1));
epi_dc_spike_lhp =mean(epi_dc_spike(BW_lhp==1));

brain_subcortica_PA = [epi_dc_stem  epi_dc_rhp  epi_dc_lhp; epi_dc_physio_stem epi_dc_physio_rhp epi_dc_physio_lhp;epi_dc_spike_stem epi_dc_spike_rhp epi_dc_spike_lhp]
