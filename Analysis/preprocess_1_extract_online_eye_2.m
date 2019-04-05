% function preprocess_1_extract_online_eye(tgt_spk_file)
% ---------
% clear
% tgt_spk_file='T017_131231_01d';
% tgt_spk_file='T011_130418_05LFP';
% tgt_spk_file='T011_130626_01b';
% tgt_spk_file='T029_150223_01b';
% tgt_spk_file='T027_140902_02b';
% tgt_spk_file='LA001_171028_s01';
% ---------

% --- load trial data from spkmatfile ---
path1=pwd;  
% addpath('S:\ShareFolder\MATLAB_work\mfile_common');   addpath_shogo; 
if str2num(tgt_spk_file(6:11))<151231,
    path2='S:\DATA_2012_upenn\units_mat_eye\';
else
    path2='S:\DATA_2015_baylor\units_mat_eye\';
end
    
if exist([path2 tgt_spk_file '.mat'],'file') & 0
    load([path2 tgt_spk_file]);
else
    
% info=parseMATFname(tgt_spk_file);
info=parseMATFname_shogo(tgt_spk_file);

% --- load video files ----
% video_folder=sprintf('%s\\%s\\%s','Y:\users\sohmae\backup_e7_shogo\data\video',info.animal,info.date);
% video_folder=sprintf('%s\\%s\\%s','Z:\Users\sohmae\backup_e7_data\video',info.animal,info.date);
% video_folder=sprintf('%s\\%s\\%s','E:\shogo\video',info.animal,info.date);
if strcmpi(tgt_spk_file(1),'l')
    video_folder=sprintf('%s\\%s\\%s','E:\Lyndsey\video',info.animal,info.date);
elseif strcmpi(tgt_spk_file(1),'c')
    video_folder=sprintf('%s\\%s\\%s','E:\Carlos\video',info.animal,info.date);
end

% video_folder=sprintf('%s\\%s\\%s','L:\shogo\video',info.animal,info.date);

if str2double(info.date)>130400 & exist([video_folder '\compressed'],'dir'), 
    video_folder=[video_folder '\compressed']; 
    ext1='_meta.mat' ;
else,
    ext1='.mat';
end

cd(video_folder), fls=dir([tgt_spk_file '_*.mat']);

for j=1:length(fls)
    
    vidfile=sprintf('%s\\%s_%s_%s_%03d%s',video_folder,info.animal,info.date,info.track,j,ext1);
    if ~exist(vidfile,'file'), continue, end
    
    load(vidfile,'metadata');

    trials.eye(j).time=metadata.eye.ts0+metadata.eye.ts_interval*[0:length(metadata.eye.trace)-1];
    trials.eye(j).trace=metadata.eye.trace;
    trials.eye(j).stimtype=lower(metadata.stim.type);
    trials.eye(j).stim=metadata.stim;
    
    switch lower(metadata.stim.type)
        case 'none'
            trials.eye(j).stimtime.st{1}=Inf;
            trials.eye(j).stimtime.en{1}=0;
            trials.eye(j).stimtime.cchan(1)=0;
        case 'puff'
            trials.eye(j).stimtime.st{1}=0;
            trials.eye(j).stimtime.en{1}=metadata.stim.totaltime;
            trials.eye(j).stimtime.cchan(1)=2;
        case 'electrical'       
            trials.eye(j).stimtime.st{1}=metadata.stim.e.delay;
            trials.eye(j).stimtime.en{1}=metadata.stim.e.traindur;
            trials.eye(j).stimtime.cchan(1)=1;
        case {'conditioning' 'electrocondition'}       
            trials.eye(j).stimtime.st{1}=0; % for CS
            trials.eye(j).stimtime.en{1}=metadata.stim.c.csdur;
            trials.eye(j).stimtime.cchan(1)=3;
            trials.eye(j).stimtime.st{2}=metadata.stim.c.isi;
            trials.eye(j).stimtime.en{2}=metadata.stim.c.usdur; % for US
            trials.eye(j).stimtime.cchan(2)=2;
        case 'optical'      
            trials.eye(j).stimtime.st{1}=metadata.stim.l.delay;
            trials.eye(j).stimtime.en{1}=metadata.stim.l.traindur;
            trials.eye(j).stimtime.cchan(1)=3;
    end
    metadata.eye=rmfield(metadata.eye,'trace');
    trials.eye(j).metadata=metadata;
end
cd(path1)

if 0
    % savefile=sprintf('%s\\%s\\%s\\%s_%s_%s.%s.mat',video_folder,info.animal,info.date,info.animal,info.date,info.track,'eye');
    % save(savefile,'trials')
    savefile=sprintf('%s\\%s_%s_%s.mat',path2,info.animal,info.date,info.track);
    save(savefile,'trials')
end
end

% setappdata(0,'trials',trials);
% OneTrialOfflineAna

