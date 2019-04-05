%%%% This is configuration file for each user. %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If Neuroblinks is launched from the root directory of the mouse, make a new directory for the session, otherwise leave that up to the user
cwd=regexp(pwd,'\\','split');
if regexp(cwd{end},'[A-Z]\d\d\d')  % Will match anything of the form LDDD, where L is single uppercase letter and DDD is a seq of 3 digits
    mkdir(datestr(now,'yymmdd'))
    cd(datestr(now,'yymmdd'))
end

% ------ Letter for mouse -----
path1=pwd;   
ind1=find(path1=='\');   
metadata.mouse=path1(ind1(end-1)+1:ind1(end)-1);

% --- camera settings ----
% Value is in microseconds and should be slightly less than interframe interval, e.g. 1/200*1e6-100 for 200 FPS
metadata.cam.init_ExposureTime = 1900;
metadata.cam.init_GainRaw = 12; % 12
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
    [9,  220,1,200, 20,0,0,0,1,1;...
     1,  220,1,200, 0,0,0,0, 1,0;...
     zeros(2,10)];
end


setappdata(0,'paramtable',paramtable);
setappdata(0,'metadata',metadata);
pause(0.3); 


metadata.stim.t.ONOFF=0;
metadata.stim.t.polarity=1;