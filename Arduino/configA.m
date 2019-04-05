%%%% This is configuration file for each user. %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% basedir='D:\shane\matlab\neuroblinks v 1.1';
basedir='D:\neuroblinks\Arduino';
% metadata.cam.recdurA=1000;

% == Will match anything of the form LDDD, where L is single uppercase letter and DDD is a seq of 3 digits ==
cwd=regexp(pwd,'\\','split');
if regexp(cwd{end},'[A-Z]\d\d\d')  
    mkdir(datestr(now,'yymmdd'))
    cd(datestr(now,'yymmdd'))
end

% ------ Letter for mouse -----
path1=pwd;   ind1=find(path1=='\');   metadata.mouse=path1(ind1(end-1)+1:ind1(end)-1);
% metadata.mouse='Sxxx';  % for Shane
% metadata.mouse='T';     % for shogo

% -- specify the location of bottomleft corner of MainWindow & AnalysisWindow  --
ghandles.pos_mainwin=[0,450];     ghandles.size_mainwin=[840 720]; 
ghandles.pos_anawin= [570 45];    ghandles.size_anawin=[1030 840]; 
ghandles.pos_oneanawin=[1020 45];    ghandles.size_oneanawin=[560 380];   
ghandles.pos_lfpwin= [570 45];    ghandles.size_lfpwin=[600 380];

% --- camera settings ----
% Value is in microseconds and should be slightly less than interframe interval, e.g. 1/200*1e6-100 for 200 FPS
metadata.cam.init_ExposureTime = 1900;
metadata.cam.init_GainRaw = 8; % 12
% NOTE: In the future this should be dynamically set based on pre and post time
% For now this variable isn't actually used by TDT version. 
metadata.cam.recdurA=1000;

% ------ Initial value of the conditioning table ----
% Search for per-mouse config file and load it if it exists, otherwise default to the paramtable below
mousedir=regexp(pwd,['[A-Z]:\\.*\\', metadata.mouse],'once','match');
condfile=fullfile(mousedir,'condparams.csv');

if exist(condfile)
	paramtable.data=csvread(condfile);
else
	paramtable.data=...
    [9,  220,1,200, 20,1,1;...
     1,  220,1,200, 0, 1,0;...
     zeros(2,7)];
end
 
comport={'COM6' 'COM6'};




% paramtable.data=zeros(10,7);
%  paramtable.data=...
%     [50, 500,1,200, 20, 1,1;...
%      50, 500,2,400, 20, 1,0;...
%      zeros(8,7)];




