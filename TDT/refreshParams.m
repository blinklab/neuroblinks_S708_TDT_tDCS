function refreshParams(hObject)

% Load objects from root app data
TDT=getappdata(0,'tdt');
metadata=getappdata(0,'metadata');
handles=guidata(hObject);
trials=getappdata(0,'trials');

val=get(handles.popupmenu_stimtype,'Value');
str=get(handles.popupmenu_stimtype,'String');

%metadata.now=datestr(now); No se si esto estaba en mi version??
metadata.stim.type=str{val};

stimmode=metadata.stim.type;

trials.savematadata=get(handles.checkbox_save_metadata,'Value');

% metadata.stim.t.ONOFF=str2double(get(handles.edit_tDCSon,'String'));
% metadata.stim.t.pulsewidth=str2double(get(handles.edit_epulsewidth,'String'));
% metadata.stim.t.traindur=str2double(get(handles.edit_etraindur,'String'));
% metadata.stim.t.amp=str2double(get(handles.edit_estimamp,'String'));
% metadata.stim.t.delay=str2double(get(handles.edit_estimdelay,'String'));
% metadata.stim.t.depth=str2double(get(handles.edit_estimdepth,'String'));

metadata.stim.e.freq=str2double(get(handles.edit_estimfreq,'String'));
metadata.stim.e.pulsewidth=str2double(get(handles.edit_epulsewidth,'String'));
metadata.stim.e.traindur=str2double(get(handles.edit_etraindur,'String'));
metadata.stim.e.amp=str2double(get(handles.edit_estimamp,'String'));
metadata.stim.e.delay=str2double(get(handles.edit_estimdelay,'String'));
metadata.stim.e.depth=str2double(get(handles.edit_estimdepth,'String'));

metadata.stim.l.freq=str2double(get(handles.edit_lstimfreq,'String'));
metadata.stim.l.pulsewidth=str2double(get(handles.edit_lpulsewidth,'String'));
metadata.stim.l.traindur=str2double(get(handles.edit_ltraindur,'String'));
metadata.stim.l.amp=str2double(get(handles.edit_lstimamp,'String'));
metadata.stim.l.delay=str2double(get(handles.edit_lstimdelay,'String'));
metadata.stim.l.depth=str2double(get(handles.edit_lstimdepth,'String'));
metadata.stim.l.ramptm=str2double(get(handles.edit_lramptm,'String'));

metadata.stim.c.puffdelay=str2double(get(handles.edit_puffdelay,'String'));
metadata.stim.c.puffdur=str2double(get(handles.edit_puffdur,'String'));

% --- conditioning -----
% Need conditional statement so we don't get an error if we're not doing conditioning so trial table hasn't been created
% if strcmpi(stimmode,'conditioning')
trialvars=readTrialTable(metadata.eye.trialnum1);
metadata.stim.c.csdur=trialvars(1);
metadata.stim.c.csnum=trialvars(2);
metadata.stim.c.isi=trialvars(3);
metadata.stim.c.usdur=trialvars(4);
% --- tDCS -----
% metadata.stim.t.ONOFF=trialvars(5);
% metadata.stim.t.polarity=trialvars(6);
% metadata.stim.t.amp=trialvars(7);
metadata.stim.t.amp=trialvars(5);
metadata.stim.t.rise_time=trialvars(6);
metadata.stim.t.fall_time=trialvars(7);
metadata.stim.t.tr_ID=trialvars(11);

% metadata.stim.c.ITI=str2double(get(handles.edit_ITI,'String'));
metadata.stim.c.ITI=trialvars(10);

metadata.stim.c.tonefreq=str2num(get(handles.edit_tone,'String'))*1000;
if length(metadata.stim.c.tonefreq)<2, metadata.stim.c.tonefreq(2)=0; end
metadata.stim.c.toneamp=str2num(get(handles.edit_toneamp,'String'));
if length(metadata.stim.c.toneamp)<2, metadata.stim.c.toneamp(2)=0; end

%=== 'auto off' codes ====%
% % Turns off the "continuous" button if we've reached the end of the table in case we want to limit the number of trials done
trialtable=getappdata(0,'trialtable');  [m,n]=size(trialtable);
if mod(metadata.eye.trialnum1,m) == 0 && metadata.eye.trialnum1 > 0
    ghandles=getappdata(0,'ghandles');
    handles=guidata(ghandles.maingui);
    if get(handles.toggle_continuous,'Value') == 1
        set(handles.toggle_continuous,'Value',0)
        set(handles.toggle_continuous,'String','Continuous: OFF')
    end
    disp(sprintf('\nEnd of Trial Table reached: Automatically disabled continuous mode\n'))
end
%=========================%
% end

puffsidestring={'ipsi' 'contra'};
metadata.stim.p.side_value=get(handles.radiobutton_contra,'Value');
metadata.stim.p.side=puffsidestring{metadata.stim.p.side_value+1};
metadata.stim.p.puffdelay=str2double(get(handles.edit_puffdelay,'String'));
metadata.stim.p.puffdur=str2double(get(handles.edit_puffdur,'String'));

switch lower(stimmode)
    case 'none'
        metadata.stim.totaltime=0;
    case 'puff'
        metadata.stim.totaltime=metadata.stim.p.puffdelay+metadata.stim.p.puffdur;
    case 'electrical'
        metadata.stim.totaltime=metadata.stim.e.traindur+metadata.stim.e.delay;
    case 'conditioning'
        % metadata.stim.totaltime=metadata.stim.c.isi+metadata.stim.c.usdur;
        metadata.stim.totaltime=metadata.stim.c.csdur+metadata.stim.c.usdur;    % So that same duration is recorded even if using two different ISIs (b/c CS dur is same)
    case {'optical','optocondition'}
        metadata.stim.totaltime=metadata.stim.l.traindur+metadata.stim.l.delay;
    case 'optoelectric'
        metadata.stim.totaltime=max(metadata.stim.e.traindur+metadata.stim.e.delay,...
            metadata.stim.l.traindur+metadata.stim.l.delay);      
    otherwise
        metadata.stim.totaltime=0;
        warning('Unknown stimulation mode set.');
end

metadata.cam.time(1)=trialvars(8);
% metadata.cam.time(1)=str2double(get(handles.edit_pretime,'String'));
metadata.cam.time(2)=metadata.stim.totaltime;
metadata.cam.time(3)=trialvars(9);

% --- saving params to memory for LFP online ana ---
tnum=metadata.cam.trialnum;
trials.params(tnum).depth=max(metadata.stim.e.depth,metadata.stim.l.depth);
trials.params(tnum).psde=metadata.stim.p.side_value;
    
setappdata(0,'trials',trials);
setappdata(0,'metadata',metadata);
