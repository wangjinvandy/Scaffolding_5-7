%%Script to get beta values for each condition for each trial in a given
%%ROI. Written by Chris McNorgan many moons ago (as far as JY knows)

%%%%%%%%%%%%%%%%%%%%%%%%%
% None of this works without Marsbar being turned on
%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize SPM and path
addpath(genpath('/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/test_folder/spm12'));
addpath('/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/test_folder/spm12/toolbox');
%addpath('/dors/gpc/JamesBooth/JBooth-Lab/BDL/fmriTools');
marsbar('on');

% What directory has all your subject folders? We assume that in each subject folder is
% a folder containing the SPM.mat file for that subject's 1st level analysis
rootDIR  = '/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/PhonReading_JW/raw_newpipeline/preprocessed';

% Where do you want the text file containing the betas to be written?
writeDIR  = rootDIR;

%How are your ROIs kept?
%1 = Everyone uses same ROI located in one location as is the case for
%standard analysis
%2= Everyone has their own ROI that is in a subject folder within a folder
%named with that ROI as is the case if paramObject was used
%3 = Everyone has their own ROI that is kept in the Analysis folder of the
%subject folders as is the case if PeakROI script is used
roistore = 2;

% If ROI store = 1 or 2, What directory has all your ROIS?
roi_file_root = [rootDIR '/IFG.oper_ROIs/'];

%what is your model directory?
modelDIR = ['/analysis/deweight'];

% make sure the scriptdir is in the path
addpath(pwd);

% set parameters
namerois = {'rIFG.oper_rhyme_vs_perc_T2_p1_k100'}; %put the name of the roi files here, omitting the '_rot.mat' part

%list all the subjects here
namesubjects={'5009' '5010' '5015' '5025' '5029' '5034' '5045' '5054' '5055' '5065' '5069' '5070' '5074' '5075' '5094' '5099' '5102' '5105' '5109' '5121' '5125' '5126' '5139' '5141' '5159' '5160' '5161' '5162' '5167' '5169' '5185' '5199' '5215' '5242' '5244' '5259' '5304'}; 


%%%%%%%%%%%%%%
%%% ALERT!!!!
%%%%%%%%%%%%%%
%are you using a .nii or _roi.mat file? If using ROI location 3, assume
%.img file instead of .nii
roi_is_image=0; %1 if using .nii file, 0 if using _roi.mat file

%name all conditions that were modeled. You must name (and get the results
%from) all conditions, whether you're interested in them or not. Edit your
%text file to show the info you're interested in. This must match what is
%listed in the SPM.xX.name variable (found in the model dir for each
%subject)
nameconditions={'Sn(1) rhyme'  'Sn(1) onset'  'Sn(1) unrel'  'Sn(1) perc'  ...
                'Sn(2) rhyme'  'Sn(2) onset'  'Sn(2) unrel'  'Sn(2) perc'  ...
                'Sn(3) rhyme'  'Sn(3) onset'  'Sn(3) unrel'  'Sn(3) perc'  ...
                'Sn(4) rhyme'  'Sn(4) onset'  'Sn(4) unrel'  'Sn(4) perc'};

%These entries correspond to the condition names in nameconditions. I'd like to use these for column headings.
friendlyconditions=	{'T1rhyme1' 'T1onset1' 'T1unrel1' 'T1perc1' ...
    'T1rhyme2' 'T1onset2' 'T1unrel2' 'T1perc2' ...
    'T2rhyme1' 'T2onset1' 'T2unrel1' 'T2perc1' ...
    'T2rhyme2' 'T2onset2' 'T2unrel2' 'T2perc2' };
    
    %If doing PPI analysis, these should look something like this
%nameconditions={'Sn(1) seed','Sn(1) Oplus_vs_Ominus','Sn(1) seedxOplus_vs_Ominus','Sn(2) seed','Sn(2) Oplus_vs_Ominus','Sn(2) seedxOplus_vs_Ominus'};
%friendlyconditions=	{'T1_seed', 'T1_contrast' 'T1_interaction' 'T2_seed', 'T2_contrast' 'T2_interaction' };

%What do you want the prefix of the text file to be? Text file will be your
%prefix+the roi name
fprefix='betas_'; %it's also possible to make our results file have a prefix.


%%%%%%%%%%%%%%%%%%should not need to edit below this line 

numsubjects = length(namesubjects);
numconditions = length(nameconditions);

    %create an roi loop
roi = 1:length(namerois);
for w = roi
    thisroi = namerois(w);
    cd(writeDIR);

    fextension='.txt'; %our results file will have a .txt extension
    writefile=char([char(fprefix) char(thisroi) char(fextension)]); %our results file will be named after the ROI associated with the betas it contains.

    delete(writefile); %deletes the data results file if it already exists (i.e., if you ran this script once already and made a whoopsie)
    headings=sprintf('%s\t',friendlyconditions{:}); % making a tab-delimited string of the condition headings in the order the betas are collected
    headings = ['ID	' headings];%We also need a column for the Subject ID
    fid=fopen(writefile, 'w');
    fprintf(fid,'%s',headings); %Now this results file has a header line containing human-readable condition names.	
    fprintf(fid,'\n');%we need a trailing newline character so our betas start on the next line. 
    % create a subject loop
    subj = 1:numsubjects;
    for x = subj
    thisguy = namesubjects(x);
    
    D=[];
    R=[];
    Y=[];
    xCON=[];
    E=[];
        fprintf('Working on participant %s for ROI %s\n', char(thisguy), char(thisroi));
	
        if roistore == 1
            if (roi_is_image)
            roi_file = char([char(roi_file_root) filesep char(thisroi) char('.nii')]);  %IF .NII FILE
            R = maroi_image(roi_file);
            else
        	roi_file = char([char(roi_file_root) filesep char(thisroi) char('_roi.mat')]); %IF .MAT FILE
            R  = maroi(roi_file);
		
            end
        elseif roistore== 2
            if (roi_is_image)
            roi_file = char([char(roi_file_root) filesep char(thisguy) filesep char(thisroi) char('.nii')]);  %IF .NII FILE
            R = maroi_image(roi_file);
            else
            roi_file = char([char(roi_file_root) filesep char(thisguy) filesep char(thisroi) char('_roi.mat')]); %IF .MAT FILE
            R  = maroi(roi_file);	
            end
        elseif roistore ==3
            roi_file = char([char(rootDIR) filesep char(thisguy) filesep modelDIR char(thisroi) char('.img')]); % ROI path
            R = maroi_image(roi_file);
        end

	

        % piece together the name of the subject directory containing the SPM.mat file
        swd = [rootDIR filesep char(namesubjects(x)) filesep modelDIR];
        %change to the subjects directory
        cd(swd);

    %try
    spm_name=load(fullfile(swd,'SPM.mat'));
    load('SPM.mat');
    cnames = transpose(SPM.xX.name); %get the condition names	

    % Make marsbar design object
    D  = mardo(spm_name);
   

    % Fetch data into marsbar data object
    Y  = get_marsy(R, D, 'mean');

    % Estimate design on ROI data
    E = estimate(D, Y);

    % get design betas
    B = betas(E);

    %catch
    %    fprintf('Encountered an error for participant %s for ROI %s\n', char(thisguy), char(thisroi));
	
%	B=zeros(1,length(cnames));
 %   end

    C=[];
    %for each condition
    %average the B value across runs
    for c = 1:(numconditions),
        betasum = 0;
        ctr = 0; %count number of betas you've grabbed for this condition. Initialize to zero because we haven't started grabbing them just yet.
        thiscondition = char(nameconditions(c));%get name of the condition as a character array
        for bindex = 1:length(B),
        %go through cnames vector
            betaname = cnames(bindex);
            found = strfind(cnames(bindex), thiscondition); %looking for whether this beta matches thiscondition
            found = length(found{1,1});
            if found > 0 %condition matches
                betasum = betasum + B(bindex); %tack on the beta for this condition, which is a Run1 beta
		ctr = ctr+1; %increment the ctr
		%If you look at the SPM.xX.name vector, you will see that the Run2 entries are 6 entries after
		%the corresponding Run1 entry. That's because there are 6 regressors per condition in the model 
		%betasum = betasum + B(bindex+6); %add on the corresponding Run2 beta, which is 6 entries further on (bindex+6)
		%ctr = ctr+1; %increment the ctr. We grabbed 2 betas, so the counter should now be at 2.                
            end
        end
    
	C(length(C)+1)=betasum/ctr;%store average value for this condition
    end

    C = [str2num(char(thisguy)) C]; %tack on thisguy's subject ID at the beginning of line
    
    cd(writeDIR);
       fmt=['%d\t', repmat('%f\t', 1, (size(C, 2)-1)), '\n'];
       fprintf(fid, fmt, transpose(C));
   % dlmwrite(writefile, C, 'delimiter', '\t', '-append');
    end
    fclose(fid);
end


