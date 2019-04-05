function TriggerStim(hObject, handles)
% this components come from pushbutton_stim_Callback of MainWindow.m

% Get stim params and pass to TDT

sendParamsToTDT(hObject)

TDT=getappdata(0,'tdt');
vidobj=getappdata(0,'vidobj');
metadata=getappdata(0,'metadata');
src=getappdata(0,'src');

metadata.TDTtankname=TDT.GetTankName();
stimmode=metadata.stim.type;
pre=metadata.cam.time(1);

if TDT.GetSysMode == 0             
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'),
    disp('%%%% TDT is Idle mode. Trigger was canceled. %%%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    return
end

% Set up camera to record
frames_per_trial=ceil(metadata.cam.fps.*(sum(metadata.cam.time))./1000);
vidobj.TriggerRepeat = frames_per_trial-1;

TDT.SetTargetVal('ustim.TrialNum',metadata.eye.trialnum2);
        
if get(handles.checkbox_record,'Value') == 1   
    % Send TDT current trial number to make mark
    TDT.SetTargetVal('ustim.CamTrial',metadata.cam.trialnum); % this will be saved in TDT storage.
    
    if get(handles.toggle_continuous,'Value') == 1,  % when continuous mode
        if TDT.GetSysMode < 3
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'),
            disp('%%%% TDT is not recording mode. Frame timings will not be saved. %%%%'),
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        end
    else
        % Make sure user knows if TDT isn't recording b/c frame times won't be recorded
        if TDT.GetSysMode < 3
            button=questdlg('You are not recording a TDT block so camera frame times will not be saved. Do you want to record a TDT block?',...
                'No active TDT block','Continue anyway','Cancel','Continue anyway');       
            switch button
                case 'Continue anyway'
                    % Do nothing and it will continue by itself
                case 'Cancel'
                    return,    % Exit stim callback
            end
        end
    end
    vidobj.StopFcn=@endOfTrial;
%     incrementStimTrial()
else
    TDT.SetTargetVal('ustim.CamTrial',0);   % Send TDT trial number of zero 
    vidobj.StopFcn=@endOfTrial;  
end

flushdata(vidobj); % Remove any data from buffer before triggering
if isprop(src,'FrameStartTriggerSource')
    src.FrameStartTriggerSource = 'FixedRate';  % Switch from free run to TTL mode
else
    src.TriggerSource = 'FixedRate';
end
start(vidobj)

metadata.ts(2)=etime(clock,datevec(metadata.ts(1)));
TDT.SetTargetVal('ustim.MatTime',metadata.ts(2));

%%%%%%%% trigger first to start camera (and second to start trial, within TDT) %%%%%%%%
TDT.SetTargetVal('ustim.StartCam',1);
pause(pre./1e3);
TDT.SetTargetVal('ustim.StartCam',0);

% if strcmpi(stimmode,'none')
%     % If doing no stim or puff
%     % Emulate button press in OpenEx
%     TDT.SetTargetVal('ustim.StartCam',1);
%     pause(pre./1e3);
%     TDT.SetTargetVal('ustim.StartCam',0);
% end
% 
% if strcmpi(stimmode,'puff')
%     TDT.SetTargetVal('ustim.PuffManual',1);
%     if get(handles.checkbox_RX6,'Value'),
%         TDT.SetTargetVal('Stim.PuffManual',1);
%         TDT.SetTargetVal('ustim.StartCam',1);
%     end
%     pause(0.01);
%     TDT.SetTargetVal('ustim.PuffManual',0);
%     if get(handles.checkbox_RX6,'Value'),
%         TDT.SetTargetVal('Stim.PuffManual',0);
%         TDT.SetTargetVal('ustim.StartCam',0);
%     end
% end
% 
% if strcmpi(stimmode,'conditioning')
%     TDT.SetTargetVal('ustim.TrigCond',1);
%     if get(handles.checkbox_RX6,'Value'),
%         TDT.SetTargetVal('Stim.TrigCond',1);
%         TDT.SetTargetVal('ustim.StartCam',1);
%     end
%     pause(0.01);
%     TDT.SetTargetVal('ustim.TrigCond',0);
%     if get(handles.checkbox_RX6,'Value'),
%         TDT.SetTargetVal('Stim.TrigCond',0);
%         TDT.SetTargetVal('ustim.StartCam',0);
%     end
% end
% 
% if strcmpi(stimmode,'optocondition')
%     TDT.SetTargetVal('ustim.StartPulse',1);
%     TDT.SetTargetVal('ustim.PuffManual',1);
% %     if get(handles.checkbox_puff,'Value')==1
%         if get(handles.checkbox_RX6,'Value'),
%             TDT.SetTargetVal('Stim.PuffManual',1);
%         end
% %     end
%     pause(0.01);
%     if get(handles.checkbox_RX6,'Value'),
%         TDT.SetTargetVal('Stim.PuffManual',0);   
%     end
%     TDT.SetTargetVal('ustim.PuffManual',0);  
%     TDT.SetTargetVal('ustim.StartPulse',0);
% end
% 
% if strcmpi(stimmode,'electrical') || strcmpi(stimmode,'optical') ||strcmpi(stimmode,'optoelectric')
%     % Emulate button press in OpenEx
%     TDT.SetTargetVal('ustim.StartPulse',1);
%     pause(0.01);
%     TDT.SetTargetVal('ustim.StartPulse',0);
% end

% --- required to initialize the eye monitor and count trial # ---- 
TDT.SetTargetVal('ustim.InitTrial',1);
pause(0.01);
TDT.SetTargetVal('ustim.InitTrial',0);

setappdata(0,'metadata',metadata);

% --- puff side swhitching ----
if strcmpi(stimmode,'puff')
    if get(handles.checkbox_puffside,'Value')
        if get(handles.radiobutton_ipsi,'Value')
            set(handles.radiobutton_contra,'Value',1)
        else
            set(handles.radiobutton_ipsi,'Value',1)
        end
    end
end


% NOTE: had to temporarily comment out this part because it causes problems with the "auto off" code that turns off
%       continuous mode when we reach the end of the trial table - SHANE
%       shogo think this issue was solved.

% ---- display current trial data in conditioning ----
if strcmpi(metadata.stim.type,'conditioning')
    
    trialvars=readTrialTable(metadata.eye.trialnum1+1);
    csdur=trialvars(1);
    csnum=trialvars(2);
    isi=trialvars(3);
    usdur=trialvars(4);
    cstone=str2num(get(handles.edit_tone,'String'));
    e_amp=str2double(get(handles.edit_estimamp,'String'));
    if length(cstone)<2, cstone(2)=0; end
    
    % --- reset background color ---
    bckgrd_color1=[1 1 1]*240/255;
    set(handles.uipanel_el,'BackgroundColor',bckgrd_color1); % light blue
    set(handles.text1,'BackgroundColor',bckgrd_color1); % light blue
    set(handles.text2,'BackgroundColor',bckgrd_color1); % light blue
    set(handles.text4,'BackgroundColor',bckgrd_color1); % light blue
    
    str2=[];   bckgrd_color2=[248 220 220]/255;
    if ismember(csnum,[5 6]), 
        str2=[' (' num2str(cstone(csnum-4)) ' Hz)'];
    elseif ismember(csnum,[7 9]), 
        str2=[' (' num2str(e_amp) ' uA)'];
        set(handles.uipanel_el,'BackgroundColor',bckgrd_color2); % light blue
        set(handles.text1,'BackgroundColor',bckgrd_color2); % light blue
        set(handles.text2,'BackgroundColor',bckgrd_color2); % light blue
        set(handles.text4,'BackgroundColor',bckgrd_color2); % light blue
    end
    
    str1=sprintf('Next:  [No %d] ID %d, ISI %d, US %d', metadata.eye.trialnum1+1, metadata.stim.t.tr_ID, isi, usdur);
    set(handles.text_disp_cond,'String',str1)
end





