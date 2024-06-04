% calculates a map of local echo times from GE field map and EPI sequence parameters
% according to analytical Eq. (1) from N. Chen et al. NeuroImage 31 (2006) 609? 622
%
% code written by Barbara Dymerska 
% code reviewed by Nadine Graedel
%
% If used please cite both:
% N. Chen et al. NeuroImage 31 (2006) 609-622
% B. Dymerska et al. Investigative Radiology 54.6 (2019): 340-348.

%%%%%%%%%%% USER PARAMETERS %%%%%%%%%%%
close all
clear
fffff
addpath('/export/home/vmalekian/Malab_lib/NIfTI_20140122/')

%% PA direction
root_dir ='/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/cmrr_mbep2d_1point25_hiedi_GRE_78_PA_0012/'; % this is also where the local TE maps will be saved
dir_GEFM= '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/distortion_rF_1/';
GEFM_name = 'topup_field.nii.gz'; % unwrapped gradient echo fieldmap in Hz
dir_epi = '/export/home/vmalekian/LC_fMRI/20240209.M700860_FIL/cmrr_mbep2d_1point25_hiedi_GRE_78_PA_0012/';
epi_name = 'PA_mc_dc_hp_physio.nii.gz';
PE_dir = +1 ; % 1 for Posterior-Anterior and -1 for Anterior-Posterior
epi_name_sample = 'PA_mc_dc_tmean.nii.gz';

% EPI parameters:
Tesp = 0.00078; % the effective echo spacing in EPI in sec
iPAT = 3 ; % total acceleration factor in phase encoding direction
pF_late = 1 ; % partial fourier last bit of k-space
pF_early = 6/8  ; % partial fourier first bit of k-space
TE = 0.019 ; % the effective echo time in seconds set as sequence parameter
my = 160 ; % matrix size in phase encoding direction

flag_voxels_out = 'yes' ; % select 'yes' or 'no', if 'yes' flags voxels with signal which will go out of k-space in EPI with value -1000

%%%%%%%%%%% END OF USER PARAMETERS %%%%%%%%%%%
output_dir = root_dir ;


disp('calculating the effective echo spacing') ;
Tesp = Tesp/iPAT ;

disp('calculating maximum/minimum TE allowed, above which type II loses occur (signal goes out of the acquired k-space)') ;
% we take into account the time it takes for excitation, navigators and any fill time. 
TE_max = TE + my*(pF_late-1/2)*Tesp;   
TE_min = TE - my*(pF_early-1/2)*Tesp; 


a1m=dir([dir_epi,epi_name]);
nf=load_untouch_nii([dir_epi,a1m(1).name]);
data= (nf.img);
EPI_spm= data;
a1m=dir([dir_GEFM,GEFM_name]);
nf=load_untouch_nii([dir_GEFM,a1m(1).name]);
GEF= (nf.img);
% GEFMrot = GEF./6.28; %% rad to Hz
GEFMrot = GEF./1; %% rad to Hz

%% VSM and Jacobian Calculation

VSM = GEFMrot*my*Tesp*PE_dir ;
% creating warp field
% VSM must be added to the PE-direction of the warp field - here 2nd dim
[Mx, My, Mz,Mt] = size(data); % reading matrix size of the EPI
[x,y,z] = ndgrid(1:Mx,1:My,1:Mz); % creating rectangular grid for 3D space
F_warp = [x(:) y(:) z(:)]; % "identity" warp field, i.e. would not do anything to our data
F_warp(:,2) = F_warp(:,2) + reshape(VSM,[Mx*My*Mz,1]); % adding VSM to war field along PE-direction

% Jacobian determinant: We have significant distortions only in
% PE-direction, so J_field boils down to the calculation of a gradient in PE-direction
VSM_grad = (circshift(VSM,1,2)-circshift(VSM,-1,2))/2 ;
VSM_grad(~isfinite(VSM_grad))=0;
J_field = 1 + VSM_grad ;
string_compcor = [dir_epi 'Jacobian.nii.gz'];
save_nii(make_nii((J_field)),string_compcor)
unix(['fslcpgeom ' dir_epi epi_name_sample ' ' string_compcor])


% disp('loading and reslicing field map to match EPI space (e.g. for oblique slices)')
% GEFM_file = fullfile(root_dir,GEFM_name) ;
% GEFM_spm = nifti(GEFM_file);
% GEFM_V = spm_vol(GEFM_file);
% data_dim = size(GEFM_spm.dat) ;   
% EPI_spm = nifti(EPI_file) ;
% 
%     GE2EPI_mat = GEFM_spm.mat\EPI_spm.mat ;
%     GEFMrot = zeros(data_dim) ;
%     data_dim_xy = data_dim(1:2);
%     
%     for slice = 1 : data_dim(3)
%         GEFMrot(:,:,slice) = spm_slice_vol(GEFM_V, GE2EPI_mat*spm_matrix([0 0 slice]), data_dim_xy, -7) ;
%     end

disp('calculating B0 gradient in EPI PE1-direction')
[G_pe1,~,~] = gradient(PE_dir*GEFMrot) ;
% Gy_file = fullfile(root_dir, 'GradPE1_GEFM.nii') ;
% createNifti(G_pe1, Gy_file, GEFM_spm.mat);

string_compcor = [dir_epi 'GradPE1_GEFM.nii.gz'];
save_nii(make_nii((G_pe1)),string_compcor)
unix(['fslcpgeom ' dir_epi epi_name_sample ' ' string_compcor])


disp('calculating local effective TE') ;
TE_mat = repmat(TE, size(G_pe1)) ;
my_mat = repmat(1/my, size(G_pe1)) ;
TEloc = TE_mat + (-G_pe1*TE*Tesp./(my_mat+G_pe1.*Tesp)) ;

if strcmp(flag_voxels_out,'yes')
    disp('flagging of the values with type II signal loss (flag = -1000)') ;
    TEloc(TEloc<TE_min) = -1000 ;
    TEloc(TEloc>TE_max) = -1000 ;

    TE_flags = zeros(size(TEloc)) ;
    TE_flags(TEloc==-1000) = 1 ;

    disp('saving flag mask') ;
%     TE_flags_file = fullfile(root_dir, 'TE_flags.nii') ;
%     createNifti(TE_flags, TE_flags_file, EPI_spm.mat);
    string_compcor = [dir_epi 'TE_flags.nii.gz'];
    save_nii(make_nii((TE_flags)),string_compcor)
    unix(['fslcpgeom ' dir_epi epi_name_sample ' ' string_compcor])

end

% TEloc_file = fullfile(root_dir, 'TE_local.nii') ;
% createNifti(TEloc, TEloc_file, EPI_spm.mat);
string_compcor = [dir_epi 'TE_local.nii.gz'];
save_nii(make_nii((TEloc)),string_compcor)
unix(['fslcpgeom ' dir_epi epi_name_sample ' ' string_compcor])

disp('saving local TE different to nominal') ;
% TEloc_file = fullfile(root_dir, 'TE_local_diff_ms.nii') ;
% createNifti((TEloc-TE)*1000, TEloc_file, EPI_spm.mat);
string_compcor = [dir_epi 'TE_local_diff_ms.nii.gz'];
save_nii(make_nii(((TEloc-TE)*1000)),string_compcor)
unix(['fslcpgeom ' dir_epi epi_name_sample ' ' string_compcor])

GEFM_name = 'jac_01.nii.gz'; % unwrapped gradient echo fieldmap in Hz
a1m=dir([dir_GEFM,GEFM_name]);
nf=load_untouch_nii([dir_GEFM,a1m(1).name]);
GEF= (nf.img);
string_compcor = [dir_epi 'jac_01_mod.nii.gz'];
save_nii(make_nii(GEF),string_compcor)
unix(['fslcpgeom ' dir_epi epi_name_sample ' ' string_compcor])






