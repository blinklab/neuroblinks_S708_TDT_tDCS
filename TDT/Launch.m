function Launch(rig,cam)
% Note that "rig" argument is not currently used for TDT option
metadata=getappdata(0,'metadata');
neuroblinks_config;	% Per user settings
configure;

setappdata(0,'paramtable',paramtable);
setappdata(0,'metadata',metadata);

%% Initialize Camera
InitCam(cam,metadata.cam.recdurA); % src and vidobj are now saved as root app data so no global vars

%% Initialize TDT interface
ghandles.TDTfig=figure;
TDT=actxcontrol('TDevAcc.X', [0 0 0 0],ghandles.TDTfig); 
set(ghandles.TDTfig,'Visible','off');     
TDT.ConnectServer('Local'); %Open a TDT Connection

ghandles.TTXfig=figure;
TTX = actxcontrol('TTank.X', [0 0 0 0],ghandles.TTXfig);
set(ghandles.TTXfig,'visible','off')


%% Set TDT Params
if TDT.GetSysMode < 1
    ok=TDT.SetTankName(tank);
    if ~ok
        disp('Cannot set tank name. Is OpenWorkbench running?');  i=0;
        while i<25 & ~ok
            disp('The system try to find TDT, again....')
            pause(1)
            
            ghandles.TDTfig=figure;
            TDT=actxcontrol('TDevAcc.X', [0 0 0 0],ghandles.TDTfig);
            set(ghandles.TDTfig,'Visible','off');
            TDT.ConnectServer('Local'); %Open a TDT Connection
            
            ghandles.TTXfig=figure;
            TTX = actxcontrol('TTank.X', [0 0 0 0],ghandles.TTXfig);
            set(ghandles.TTXfig,'visible','off')
            
            ok=TDT.SetTankName(tank);
            i=i+1;
        end
        if ~ok
            error('Cannot set tank name. Is OpenWorkbench running?')
        end
    end
    ok=TDT.SetSysMode(2);
    if ~ok
    	error('Cannot set TDT to Preview mode. Is OpenWorkbench running?')
    end
else 
    tank=TDT.GetTankName();
end

if TTX.ConnectServer('Local','Me') == 0, error('Error connecting to server'); end
if TTX.OpenTank(tank,'R') == 0, error('Error opening tank'); end


%% Save objects to root app data
setappdata(0,'tdt',TDT)
setappdata(0,'ttx',TTX)

%% Set up timer for eyelid streaming

% 200 Hz timer
% First delete old instances
% t=timerfind('Name','eyelidTimer');
% delete(t)
% eyelidTimer=timer('Name','eyelidTimer','Period',0.005,'ExecutionMode','FixedRate','TimerFcn',@eyelidstream,'BusyMode','drop');

%% Open GUI
ghandles.maingui=MainWindow;
% movegui(ghandles.maingui,pos_mainwin)
set(ghandles.maingui,'units','pixels')
ghandles.pos_mainwin = ghandles.pos_mainwin+(rig-1)*ghandles.shift_for_2nd_rig;
ghandles.pos_oneanawin = ghandles.pos_oneanawin+(rig-1)*ghandles.shift_for_2nd_rig;
set(ghandles.maingui,'position',[ghandles.pos_mainwin ghandles.size_mainwin])

% Save handles to root app data
setappdata(0,'ghandles',ghandles)



