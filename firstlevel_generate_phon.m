% First level analysis, wrote by Jin Wang 3/15/2019
% You should define your conditions, onsets, duration, TR.
% The repaired images will be deweighted from 1 to 0.01 in the first level
% estimation (to changed the default art_redo.m using art_deweight.txt to using art_repaired.txt).
% The 6 movement parameters we got from realignment is added into the model regressors to remove the small motion effects on data.
% The spm12 we used has the following changes
% spm_defaults.m defaults.mask.thresh = 0.5 (insetead of the default value of 0.8); 
% toolbox included marsbar and art-repair;
% we modified the art_pair toolbox  art_redo.m as art_redo_jin.m to fit our code structure


clear all;
addpath(genpath('/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/PhonReading_JW/raw_newpipeline/scripts'));
spm_path='/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/test_folder/spm12';
addpath(genpath(spm_path));

%define your data path
data=struct();
root='/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/PhonReading_JW/raw_newpipeline';
subjects={'5259' '5304'};
% '5009' '5010' '5015' '5025' '5029' '5034' '5045' '5054' '5055' '5065' '5069'
% '5070' '5074' '5075' '5094' '5099' '5102' '5105' '5109' '5121' '5125' '5126' '5139'
% '5141' '5159' '5160' '5161' '5162' '5167' '5169' '5199' '5215' '5242' '5244'
% 
global CCN
CCN.session='ses*';
CCN.func_pattern='T*';
CCN.file='wtask*bold.nii';
CCN.files='wtask*bold_0*';
analysis_folder='mvpa_glm';
model_deweight='deweight';
CCN.preprocessed='preprocessed';
CCN.rpfile='rp_*.txt';

%define your task conditions, onsets
%Phon task, each run should have its own condition defined to be run.
conditions{1}={'rhyme' 'onset' 'unrelated' 'perc'};
conditions{2}={'rhyme' 'onset' 'unrelated' 'perc'};
conditions{3}={'rhyme' 'onset' 'unrelated' 'perc'};
conditions{4}={'rhyme' 'onset' 'unrelated' 'perc'};

%load onsets files, be aware of the sequence, it should be consistent with
%your conditions
a=load('/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/onsets/phon_1_rhyme.txt');
b=load('/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/onsets/phon_1_onset.txt');
c=load('/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/onsets/phon_1_unrel.txt');
d=load('/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/onsets/phon_1_perc.txt');
e=load('/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/onsets/phon_2_rhyme.txt');
f=load('/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/onsets/phon_2_onset.txt');
g=load('/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/onsets/phon_2_unrel.txt');
h=load('/dors/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/onsets/phon_2_perc.txt');
onsets{1}(:,1)=a;
onsets{1}(:,2)=b;
onsets{1}(:,3)=c;
onsets{1}(:,4)=d;
onsets{2}(:,1)=e;
onsets{2}(:,2)=f;
onsets{2}(:,3)=g;
onsets{2}(:,4)=h;
onsets{3}(:,1)=a;
onsets{3}(:,2)=b;
onsets{3}(:,3)=c;
onsets{3}(:,4)=d;
onsets{4}(:,1)=e;
onsets{4}(:,2)=f;
onsets{4}(:,3)=g;
onsets{4}(:,4)=h;

%duration
dur=0;

%TR
TR=1.25;

%define your contrasts, make sure your contrasts and your weights should be
%matched.
contrasts={'onset_vs_perc_T1' ...
    'rhyme_vs_perc_T1' ...
    'onset_vs_perc_T2' ...
    'rhyme_vs_perc_T2'};

onset_vs_perc=[0 1 0 -1];
rhyme_vs_perc=[1 0 0 -1];
empty=[0 0 0 0];

%adjust the contrast by adding six 0s into the end of each session
rp_w=zeros(1,6);
weights={[onset_vs_perc rp_w onset_vs_perc rp_w empty rp_w empty rp_w] ...
    [rhyme_vs_perc rp_w rhyme_vs_perc rp_w empty rp_w empty rp_w] ...
    [empty rp_w empty rp_w onset_vs_perc rp_w onset_vs_perc rp_w] ...
    [empty rp_w empty rp_w rhyme_vs_perc rp_w rhyme_vs_perc rp_w] };

%%%%%%%%%%%%%%%%%%%%%%%%Do not edit below here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if you define your contrasts in a correct way
if length(weights)~=length(contrasts)
    error('the contrasts and the weights are not matched');
end

% Initialize
%addpath(spm_path);
spm('defaults','fmri');
spm_jobman('initcfg');
spm_figure('Create','Graphics','Graphics');

% Dependency and sanity checks
if verLessThan('matlab','R2013a')
    error('Matlab version is %s but R2013a or higher is required',version)
end

req_spm_ver = 'SPM12 (6225)';
spm_ver = spm('version');
if ~strcmp( spm_ver,req_spm_ver )
    error('SPM version is %s but %s is required',spm_ver,req_spm_ver)
end

%Start to analyze the data from here
for i=1:length(subjects)
    fprintf('work on subject %s', subjects{i});
    CCN.subject=[root '/' CCN.preprocessed '/' subjects{i}];
    %specify the outpath,create one if it does not exist
    out_path=[CCN.subject '/' analysis_folder];
    if ~exist(out_path)
        mkdir(out_path)
    end
    
    %specify the deweighting spm folder, create one if it does not exist
    model_deweight_path=[out_path '/' model_deweight];
    if exist(model_deweight_path,'dir')~=7
        mkdir(model_deweight_path)
    end
    
    %find folders in func
    CCN.functional_dirs='[subject]/[session]/func/[func_pattern]/';
    functional_dirs=expand_path(CCN.functional_dirs);
    
    %load 6 movement parameters
    mv=[];
    rp_file=expand_path([CCN.functional_dirs '[rpfile]']);
    for i=1:length(rp_file)
        rp=load(rp_file{i});
        mv{i}=rp;
    end
    data.mv=mv;
    
    %load the functional data
    swfunc=[];
    P=[];
    for j=1:length(functional_dirs)
        P{j}=char(expand_path([functional_dirs{j} '[files]']));
        if isempty(P{j}) % if there is no expanded data, then expand the preprocessed data
            file=expand_path([functional_dirs{j} '[file]']);
            hdr = load_nii_hdr(char(file));
            nscan = hdr.dime.dim(5);
            expand_nii_scan(char(file),1:nscan);
            P{j}=char(expand_path([functional_dirs{j} '[files]']));
        end
        for ii=1:size(P{j},1)
            swfunc{j}(ii,:)=[P{j}(ii,:) ',1'];
        end
    end
    data.swfunc=swfunc;
    
    %pass the experimental design information to data
    data.conditions=conditions;
    data.onsets=onsets;
    data.dur=dur;
    
    %run the firstlevel modeling and estimation (with deweighting)
    mat=firstlevel(data, out_path, TR, model_deweight_path);
    origmat=[out_path '/SPM.mat'];
    %run the contrasts
    contrast(origmat,contrasts,weights);
    contrast(mat,contrasts,weights);
    
end