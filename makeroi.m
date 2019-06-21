%
% I should get into the habit of putting my name in these things...
% makeroi.m by Chris McNorgan on August 30, 2012
% version 1.1
%
% function makeroi(paramObject)
% 1) creates a subdirectory in the current working directory 
% 2) for each subject {s1 ... sn} creates a subject sub-subdirectory
% 3) containing ROI directories for each region {r1 ... rn}
% paramObject is a structure variable with a number of required and optional fields:
%	REQUIRED FIELDS:
%	paramObject.data_root='/path/to/dir'; %the parent directory for all the subject folders
%	paramObject.SPM_mat_folder='subdirectory/names'; %the name of the directory under the subject data folder containing the SPM.mat file with contrast info
%	paramObject.contrast_name='name_of_contrast'; %the name of the contrast associated with the statistical map you want to use. Check your SPM.mat file for names
%	paramObject.regions={'R1' ... 'Rn'}; %for each region 'Rx', there should be a tab-delimited list of xyz coordinates called Rx.txt in the working directory.
%	paramObject.images={'I1' ... 'In'}; for each filename 'In', there should be a binarized nifti file image In.nii in your current working directory.
%	paramObject.subjects={'S1' ... 'Sn'}; %for each subject 'Sx', there should be a data folder called Sx in the data_root

%	OPTIONAL FIELDS:
%	If these values are not specified, defaults will be used (p=.05, radius=5mm, k=100%)
%	paramObject.p = p; %uncorrected p-value for statistical map with which sphere is intersected
%	paramObject.radius = r; %sphere radius, in mm
%	paramObject.k = 'k'; %specify paramObject.k as either 'k%' or just 'k'. 
%			NOTE that 'k' is a STRING (character array), not a number
%			if k is specified as 'k%' (e.g., '50%'), it will take the top k percentile of voxels
%			if k is specified as k (e.g., '50') , then it will take the top k voxels
%			It will always select at least 1 voxel (e.g., top 10% of a 9 voxel blob < 1, so select the top 1 voxel)
%	paramObject.savesphere=[1/0]; %do you want to save the base spheres for each person? Might be useful for debugging failed ROI attempts
%	paramObject.savedir='path/to/dir';%name of the subdirectory to be created in the current directory into which the ROIs will be saved
%			if a savedir is not specified, a directory with today's date will be created


function  makeroi(paramObject)

	if(length(which('marsbar.m')) < 1)
	error('This function relies on the marsbar toolbox which was not found in your path');
	return
	end

	marsbar('on');
	def_p=.05;
	def_radius=5;
	def_k=1;
	conname=paramObject.contrast_name;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Set the optional parameters
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	%
	%Sphere radius
	try
		sphere_radius=double(int8(paramObject.radius));
		if(sphere_radius > 12)
			fprintf(1,'Not that I am judging but, your sphere radius of %d seems pretty big. A data type conversion problem perhaps?\n', sphere_radius); 
			
		end
	catch err
		sphere_radius=def_radius;
	end

	%
	%p-threshold for activation map	
	try
		p_thresh=double(paramObject.p);
	catch err
		fprintf(1,'No p-value threshold provided for activation map. Using default uncorrected p<%f\n', def_p);
		p_thresh=def_p;

	end
	
	%param_details: used for naming the ROIs created for spheres
	%nii_param_details: used for naming ROIS created using nii masks

	param_details=[conname '_' num2str(sphere_radius) 'mm_p' num2str(p_thresh)];
	nii_param_details=[conname '_p' num2str(p_thresh)];
	try
		param_details=[param_details '_k' strrep(paramObject.k,'%','pct')];
		nii_param_details=[nii_param_details '_k' strrep(paramObject.k,'%','pct')];
	catch err
	end

	try
		roi_dir=paramObject.savedir;
	catch err
		roi_dir=['ROI_' date()];
	end
	

	
	root_dir=pwd;
	roi_root=[root_dir filesep roi_dir];
	nsubjects=length(paramObject.subjects);

	%this next part pieces together a bunch of parameters to figure out the path to the SPM.mat files
	data_dirs=[repmat(paramObject.data_root,nsubjects,1) ...
        repmat(filesep,nsubjects,1) ...
	    char(paramObject.subjects) ...
        repmat(filesep,nsubjects,1) ...
	    repmat(paramObject.SPM_mat_folder,nsubjects,1)];
	
	if(~exist(roi_root))
		mkdir(fullfile(roi_root));
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%execute loop for each region
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for w = 1:length(paramObject.regions)
		
		brainregion = char(paramObject.regions(w));
		details=[brainregion param_details];
		centerfile = [root_dir filesep brainregion '.txt'];
    		sphere_centres = load(centerfile);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%A little complimentary error checking: there should be exactly 1 
		%triplet of coordinates per subject folder
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if(size(sphere_centres,1) ~= nsubjects)
			
			error(['ROI coordinates for ' brainregion ' mismatch number of subjects']);
			return;
		end


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%execute nested loop for each subject
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		for x=1:nsubjects

			%First, we need to generate a spherical ROI for the subject
			roi_dir=[roi_root filesep char(paramObject.subjects(x))];
			if(~exist(roi_dir))
				mkdir(fullfile(roi_dir));
			end
			sphere_centre = sphere_centres(x,:);

		        sphere_roi = maroi_sphere(struct('centre', sphere_centre, 'radius', sphere_radius));
		

		        % optional: save spherical ROI to MarsBaR ROI file, in subject's ROI
			try
				if(paramObject.savesphere == 1)
					details = [brainregion];
			       	 	roitosave = label(sphere_roi, details);
			       	 	detailsmat = [details, '_sphere_roi.mat'];
					detailsnii = [details, '_sphere_roi.nii'];
					saveroi(roitosave, fullfile(roi_dir,detailsmat ));
					save_as_image(roitosave, fullfile(roi_dir, detailsnii));
				end
			catch exception
				%do nothing
			end

			%now we have the sphere, we need a statistical map
			
 			% Get SPM model
                try
		        model_dir = data_dirs(x,:);
		        smodel = mardo(fullfile([model_dir filesep 'SPM.mat']));
                catch err
                    %The script inexplicably started inserting spaces in
                    %the directory names after the subject number. This next line strips them. In
                    %the even that your directory is SUPPOSED to have spaces in the
                    %names, you'll have to do some different hack.
                    model_dir(model_dir==' ')='';
                end
		        if ~is_spm_estimated(smodel)
          			error(['Session model has not been estimated. ' ...
            				'You may need to run the run_preprocess script']);
        		end

		        % Get activation cluster by loading T image
        		con_name = paramObject.contrast_name;
        		t_con = get_contrast_by_name(smodel, con_name);
        		if isempty(t_con)
          			error(['Cannot find the contrast ' con_name ...
            			' in the design; has it been estimated?']);
       			end
		        % SPM2 stores contrasts as vols, SPM99 as filenames
		        if isstruct(t_con.Vspm)
		        	t_con_fname = t_con.Vspm.fname;
        		else
          			t_con_fname = t_con.Vspm;
        		end

        		t_img_name = fullfile([model_dir filesep t_con_fname]);
        		if ~exist(t_img_name)
          			error(['Cannot find t image ' t_img_name ...
             				'; has it been estimated or deleted?']);
        		end

			try
			        y = getdata(sphere_roi, t_img_name); %get all voxels in the sphere (regardless of activity)
			catch err
				%in my testing, the getdata command would fail depending on the coordinates I used. My guess
				%is that you can't get coordinates outside of the brain (maybe including CSF)
				fprintf(1,'Unable to find data at %d %d %d. Are your coordinates reasonable?\n', ...
	sphere_centre(1), sphere_centre(2), sphere_centre(3));
			end
        		n_sphere_voxels=size(y,2); %how many voxels are in the sphere?

			erdf = error_df(smodel); %retrieve the error degrees of freedom (to calculate the t-threshold)
		        t_thresh = spm_invTcdf(1-p_thresh, erdf); %calculate the t-threshold for the p-value we're using

			%find the voxels in contrast image that are at the calculated t_threshold
			V = spm_vol(t_img_name);
        		img = spm_read_vols(V);
      			tmp = find(img(:) >= t_thresh);
        		img = img(tmp);
        		XYZ = mars_utils('e2xyz', tmp, V.dim(1:3));

			act_roi = maroi_pointlist(struct('XYZ', XYZ, 'mat', V.mat), 'vox');

			%simplest case: the intersection of the sphere and the activation map:	
		        masked_roi = act_roi & sphere_roi;
                	final_roi = masked_roi;
			%%%%%%%%%%%%%%%%%%%%%%%%
			%Works to this point.
			%%%%%%%%%%%%%%%%%%%%%%%

			%
			%What if we want exactly k or k% voxels in the roi?
			%
    		try
                [Y, multv, vXYZ, mat] = getdata(masked_roi, t_img_name);
                tmp=sort(Y,'descend');
                n_masked_voxels=length(tmp);
            catch err
               tmp=[Inf]
               n_masked_voxels=0;
               fprintf(1,'The %s ROI for %s is empty. Try a more liberal p-value\n', brainregion, char(paramObject.subjects(x)));
            end

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% calculate k number of voxels (optional parameter) 
			%(interpreted as a percent or number of voxels)
       			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			try
				%default: k is a fixed number
				k=str2num(paramObject.k);
				if(~isempty(findstr(paramObject.k,'%')))
				%but if k is specified with a percent sign, calculate what number that works out to be...

					%e.g., k='50%' and the masked roi has 128 voxels. So we want the top (0.50 * 128) = 64 voxels
					k=(str2num(strrep(paramObject.k,'%','')));%strip out the percent.
					k=k/100; %convert to a fraction
					k=double(int8(n_masked_voxels*k));%multiply the fraction by the number of voxels to calculate k
					if(k<1)
						k=1;
					end
				end
			
			catch err
				fprintf(1,'No voxel count specified. Using 100 percent of intersection of sphere and activation map\n');
				k=n_masked_voxels;
		
			end

			%filter to include only the top k or k% of the masked sphere
			%only if k is a subset
			if(k < n_masked_voxels)
				
				k_cut=tmp(k);
				filter=ones(size(Y));%generate a transparent mask
				filter(find(Y < k_cut))=0;%zero out the sub-threshold spots
				filtered_coords=vXYZ(:,find(filter > 0));%apply filter to the coordinates from the mask
				filtered_roi = maroi_pointlist(struct('XYZ', filtered_coords, 'mat', V.mat), 'vox'); %convert coordinates into an ROI
                		final_roi=masked_roi & filtered_roi;
			end

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% SAVE THE ROI object
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			details_nii = [brainregion '_' param_details '_roi.nii'];
			details_mat = [brainregion '_' param_details '_roi.mat'];
			roitosave = label(final_roi, [brainregion '_' param_details]);
        		% save ROI to MarsBaR ROI file, in ROI directory
			%first check if destination directory for subject exists
			if (~exist(roi_dir))
				mkdir(roi_dir); %make it if necessary
			end
        	
        		saveroi(roitosave, fullfile(roi_dir,details_mat ));
		        save_as_image(roitosave, fullfile(roi_dir, details_nii));
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% finish nested subject loop
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		end
	
	%%%%%%%%%%%%%%%%%%%%%
	% finish region loop
	%%%%%%%%%%%%%%%%%%%%	
	end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%execute loop for each nii image file
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for w = 1:length(paramObject.images)
		
		brainregion = char(paramObject.images(w));
		details=[brainregion param_details];
		imagefile = [root_dir filesep brainregion '.nii'];
    		

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%execute nested loop for each subject
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		for x=1:nsubjects

			%First, we need to generate a spherical ROI for the subject
			roi_dir=[roi_root filesep char(paramObject.subjects(x))];
			if(~exist(roi_dir))
				mkdir(fullfile(roi_dir));
            end
                
		        image_roi = maroi_image(imagefile);
		

			%now we have the ROI image, we need a statistical map
			
 			% Get SPM model

		        model_dir = data_dirs(x,:);
		        smodel = mardo(fullfile([model_dir filesep 'SPM.mat']));
		        if ~is_spm_estimated(smodel)
          			error(['Session model has not been estimated. ' ...
            				'You may need to run the run_preprocess script']);
        		end

		        % Get activation cluster by loading T image
        		con_name = paramObject.contrast_name;
        		t_con = get_contrast_by_name(smodel, con_name);
        		if isempty(t_con)
          			error(['Cannot find the contrast ' con_name ...
            			' in the design; has it been estimated?']);
       			end
		        % SPM2 stores contrasts as vols, SPM99 as filenames
		        if isstruct(t_con.Vspm)
		        	t_con_fname = t_con.Vspm.fname;
        		else
          			t_con_fname = t_con.Vspm;
        		end

        		t_img_name = fullfile([model_dir filesep t_con_fname]);
        		if ~exist(t_img_name)
          			error(['Cannot find t image ' t_img_name ...
             				'; has it been estimated or deleted?']);
        		end

			try
			        y = getdata(image_roi, t_img_name); %get all voxels in the image (regardless of activity)
			catch err
				%in my testing, the getdata command would fail depending on the coordinates I used. My guess
				%is that you can't get coordinates outside of the brain (maybe including CSF)
				fprintf(1,'Unable to find data in %s. Are your coordinates reasonable?\n', imagefile);
			end
        		n_voxels=size(y,2); %how many voxels are in the image?

			erdf = error_df(smodel); %retrieve the error degrees of freedom (to calculate the t-threshold)
		        t_thresh = spm_invTcdf(1-p_thresh, erdf); %calculate the t-threshold for the p-value we're using

			%find the voxels in contrast image that are at the calculated t_threshold
			V = spm_vol(t_img_name);
        		img = spm_read_vols(V);
      			tmp = find(img(:) >= t_thresh);
        		img = img(tmp);
        		XYZ = mars_utils('e2xyz', tmp, V.dim(1:3));

			act_roi = maroi_pointlist(struct('XYZ', XYZ, 'mat', V.mat), 'vox');

			%simplest case: the intersection of the sphere and the activation map:	
		        masked_roi = act_roi & image_roi;
                	final_roi = masked_roi;
			%%%%%%%%%%%%%%%%%%%%%%%%
			%Works to this point.
			%%%%%%%%%%%%%%%%%%%%%%%

			%
			%What if we want exactly k or k% voxels in the roi?
			%
            %this is in a try/catch block because its possible that an S
            %will not have voxels within the mask in t_img_name that are
            %above your p-value threshold. If so, handle gracefully
			try
                [Y, multv, vXYZ, mat] = getdata(masked_roi, t_img_name);
                tmp=sort(Y,'descend');
                n_masked_voxels=length(tmp);
            catch err
               tmp=[Inf]
               n_masked_voxels=0;
               fprintf(1,'The %s ROI for %s is empty. Try a more liberal p-value\n', brainregion, char(paramObject.subjects(x)));
            end
            
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% calculate k number of voxels (optional parameter) 
			%(interpreted as a percent or number of voxels)
       			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			try
				%default: k is a fixed number
				k=str2num(paramObject.k);
				if(~isempty(findstr(paramObject.k,'%')))
				%but if k is specified with a percent sign, calculate what number that works out to be...

					%e.g., k='50%' and the masked roi has 128 voxels. So we want the top (0.50 * 128) = 64 voxels
					k=(str2num(strrep(paramObject.k,'%','')));%strip out the percent.
					k=k/100; %convert to a fraction
					k=double(int8(n_masked_voxels*k));%multiply the fraction by the number of voxels to calculate k
					if(k<1)
						k=1;
					end
				end
			
			catch err
				fprintf(1,'No voxel count specified. Using 100 percent of intersection of sphere and activation map\n');
				k=n_masked_voxels;
		
			end

			%filter to include only the top k or k% of the masked image
			%only if k is a subset
			if(k < n_masked_voxels)
				k_cut=tmp(k);
				filter=ones(size(Y));%generate a transparent mask
				filter(find(Y < k_cut))=0;%zero out the sub-threshold spots
				filtered_coords=vXYZ(:,find(filter > 0));%apply filter to the coordinates from the mask
				filtered_roi = maroi_pointlist(struct('XYZ', filtered_coords, 'mat', V.mat), 'vox'); %convert coordinates into an ROI
                		final_roi=masked_roi & filtered_roi;
			end

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% SAVE THE ROI object
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			details_nii = [brainregion '_' nii_param_details '_roi.nii'];
			details_mat = [brainregion '_' nii_param_details '_roi.mat'];
			roitosave = label(final_roi, [brainregion '_' nii_param_details]);
        		% save ROI to MarsBaR ROI file, in ROI directory
			%first check if destination directory for subject exists
			if (~exist(roi_dir))
				mkdir(roi_dir); %make it if necessary
			end
        	
        		saveroi(roitosave, fullfile(roi_dir,details_mat ));
		        save_as_image(roitosave, fullfile(roi_dir, details_nii));
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% finish nested subject loop
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		end
	
	%%%%%%%%%%%%%%%%%%%%%
	% finish region loop
	%%%%%%%%%%%%%%%%%%%%	
	end

end

