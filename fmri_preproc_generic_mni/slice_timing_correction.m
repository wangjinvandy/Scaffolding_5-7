function afunc_file = slice_timing_correction(func_file,params)
% Slice timing correction is applied in the third dimension (Z, often
% axial).


% Expand func filenames. A single 4D Nifti file is expected per run. An
% SPM-style list of volumes is returned.
[func_p,func_n,func_e] = fileparts(func_file);
func_vols = cellstr(spm_select('ExtFPList',func_p,['^' func_n func_e '$'],inf));


% Number of slices
V = spm_vol(func_vols{1});
nslices = V.dim(3);


% Slice order
switch params.slorder
	
	case 'ascending'
		slorder = 1:nslices;
		
	case 'descending'
		slorder = nslices:-1:1;

	case 'ascending_interleaved'
		slorder = [1:2:nslices 2:2:nslices];
		
	case 'descending_interleaved'
		slorder = [nslices:-2:1 (nslices-1):-2:1];

	case {'none',''}
		warning('Skipping slice timing correction')
		[p,n,e] = fileparts(func_file);
		afunc_file = fullfile(p,['a' n e]);
		copyfile(func_file,afunc_file);
		return
		
	otherwise
		error('Unknown slice order')
		
end


% Batch job for slice timing correction
matlabbatch = [];
tag = 1;
matlabbatch{tag}.spm.temporal.st.scans = {func_vols};
matlabbatch{tag}.spm.temporal.st.nslices = nslices;
matlabbatch{tag}.spm.temporal.st.tr = params.tr;
matlabbatch{tag}.spm.temporal.st.ta = params.tr - params.tr/nslices;
matlabbatch{tag}.spm.temporal.st.so = slorder;
matlabbatch{tag}.spm.temporal.st.refslice = 1;
%save(fullfile(func_p,'batch_slice_timing_correction.mat'),'matlabbatch')
spm_jobman('run',matlabbatch)


% Filename for slice time corrected images. The file doesn't
% exist, so we have to predict its name
clear afunc_filename
[func_p,func_n,func_e] = fileparts(func_file);
afunc_file = fullfile(func_p,['a' func_n func_e]);

