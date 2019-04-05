function sendParamsToTDT(hObject)

refreshParams(hObject)

% Load objects from root app data
TDT=getappdata(0,'tdt'); 
metadata=getappdata(0,'metadata');

handles=guidata(hObject);

% Pass pulse values for Camera to TDT - even if we're not actually recording to disk
T_length=(sum(metadata.cam.time))./1000;
frames_per_trial=ceil(metadata.cam.fps.*T_length);
TDT.SetTargetVal('ustim.FramePulse',1e3/(2*metadata.cam.fps));
TDT.SetTargetVal('ustim.NumFrames',frames_per_trial);

% ---- for electrical stimulation ----
TDT.SetTargetVal('ustim.EPulseFreq',metadata.stim.e.freq);
TDT.SetTargetVal('ustim.EPulseWidth',metadata.stim.e.pulsewidth);
TDT.SetTargetVal('ustim.ETrainDur',metadata.stim.e.traindur);
TDT.SetTargetVal('ustim.EStimAmp',metadata.stim.e.amp);
TDT.SetTargetVal('ustim.EStimDelay',metadata.stim.e.delay);
TDT.SetTargetVal('ustim.EStimDepth',metadata.stim.e.depth);

% ---- for Laser stimulation ----
TDT.SetTargetVal('ustim.LPulseFreq',metadata.stim.l.freq);
TDT.SetTargetVal('ustim.LPulseWidth',metadata.stim.l.pulsewidth);
TDT.SetTargetVal('ustim.LTrainDur',metadata.stim.l.traindur);
TDT.SetTargetVal('ustim.LStimAmp',metadata.stim.l.amp);
TDT.SetTargetVal('ustim.LStimDelay',metadata.stim.l.delay);
TDT.SetTargetVal('ustim.LStimDepth',metadata.stim.l.depth);
TDT.SetTargetVal('ustim.RampTime',metadata.stim.l.ramptm);

% ---- for tDCS stimulation ----
TDT.SetTargetVal('ustim.tDCS_rise_T',metadata.stim.t.rise_time); %tDCS on=1/off=0
TDT.SetTargetVal('ustim.tDCS_fall_T',metadata.stim.t.fall_time); %tDCS on=1/off=0
TDT.SetTargetVal('ustim.TDCSStimAmp',metadata.stim.t.amp); %tDCS on=1/off=0
% metadata.stim.t.tr_ID
% TDT.SetTargetVal('ustim.TDCSPulseDur',metadata.stim.t.pulsewidth);
% TDT.SetTargetVal('ustim.Polarity',metadata.stim.t.traindur); %anodal=0/catodal=1
% TDT.SetTargetVal('ustim.TDCSStimAmp',metadata.stim.t.amp);
% TDT.SetTargetVal('ustim.TDCSPauseDur',metadata.stim.t.delay); %Tiempo desde que empieza el tDCS hasta que estimulamos
% TDT.SetTargetVal('ustim.LStimDepth',metadata.stim.t.depth);
% TDT.SetTargetVal('ustim.TDCSRampTime',metadata.stim.t.ramptm);

if metadata.stim.l.ramptm > 0
    TDT.SetTargetVal('ustim.PulseShape',1);
else
    TDT.SetTargetVal('ustim.PulseShape',0);
end

% Have to subtract 1 b/c TDT is zero-referenced and Matlab is
% one-referenced

% --- conditioning -----
% Need conditional statement so we don't get an error if we're not doing conditioning so trial table hasn't been created
% if strcmpi(metadata.stim.type,'conditioning')
    % --- behavioral stim by RZ5 ---- 
    TDT.SetTargetVal('ustim.ITI',metadata.stim.c.ITI);
    TDT.SetTargetVal('ustim.CsDur',metadata.stim.c.csdur);
    csnum=metadata.stim.c.csnum;  cstonefreq=0;  cstoneamp=0;  csnum1=3; csnum2=3;
    if ismember(csnum,[5 6]),  % for auditory CS
        cstonefreq=min(metadata.stim.c.tonefreq(csnum-4), 40000);  
        cstoneamp=metadata.stim.c.toneamp(csnum-4);
        csnum=0; csnum1=0; csnum2=3; 
    elseif ismember(csnum,[1 2 3 4]),  % for DIO CS (LED/Wisker)
        csnum1=1; csnum2=csnum-1;
    end
    if strcmpi(metadata.stim.type,'conditioning') & ismember(csnum,[7 9]),  % for electrical CS
        TDT.SetTargetVal('ustim.ETrainDur',metadata.stim.c.csdur);
        TDT.SetTargetVal('ustim.EStimDelay',0);
        csnum1=3; csnum2=3;
    end
    if strcmpi(metadata.stim.type,'conditioning') & ismember(csnum,[8 9]),  % for Laser CS
        TDT.SetTargetVal('ustim.LTrainDur',metadata.stim.c.csdur);
        if metadata.stim.l.freq<=2,
            TDT.SetTargetVal('ustim.LPulseWidth',metadata.stim.c.csdur*1000);
        end
        TDT.SetTargetVal('ustim.LStimDelay',0);
        csnum1=3; csnum2=3;
    end
    TDT.SetTargetVal('ustim.CSNum',csnum);
    TDT.SetTargetVal('ustim.CSNum1',csnum1); % 0 tone, 1 DIO output (LED/wisker), 3 el/opt
    TDT.SetTargetVal('ustim.CSNum2',csnum2); % for DIO, 0 DIO1, 1 DIO2, 2 DIO3 
    TDT.SetTargetVal('ustim.CSToneFreq',cstonefreq);
    TDT.SetTargetVal('ustim.CSToneAmp',cstoneamp);
    TDT.SetTargetVal('sound.CSToneFreq',cstonefreq);
    TDT.SetTargetVal('sound.CSToneAmp',cstoneamp);
    % csnum,cstonefreq,

    TDT.SetTargetVal('ustim.ISI',metadata.stim.c.isi);
    TDT.SetTargetVal('ustim.UsDur',metadata.stim.c.usdur);
    % TDT.SetTargetVal('ustim.PuffMDelay',metadata.stim.c.puffdelay);
    TDT.SetTargetVal('ustim.PuffDurM',metadata.stim.c.puffdur);
% end

TDT.SetTargetVal('ustim.PreStimTime',metadata.cam.time(1));
TDT.SetTargetVal('ustim.PuffSide',metadata.stim.p.side_value);
TDT.SetTargetVal('ustim.PuffMDelay',metadata.stim.p.puffdelay);


% switch lower(metadata.stim.type)
%     case 'electrical'
%         TDT.SetTargetVal('ustim.StimDevice',0);
%     case {'optical','optocondition'}
%         TDT.SetTargetVal('ustim.StimDevice',1);
%     case 'optoelectric'
%         TDT.SetTargetVal('ustim.StimDevice',2);
% end

switch lower(metadata.stim.type)
    case {'conditioning','electrocondition','optocondition'}
        TDT.SetTargetVal('ustim.TrialType',0);
    case {'puff'}
        TDT.SetTargetVal('ustim.TrialType',1);
    case {'none','electrical','optical','optoelectric'} 
        TDT.SetTargetVal('ustim.TrialType',3);   % no output 
end

switch lower(metadata.stim.type) % controling el or opt devices
    case {'none','puff'}
        TDT.SetTargetVal('ustim.StimDevice',3);   % no output from el/opt devices
    case {'electrical'}
        TDT.SetTargetVal('ustim.StimDevice',0);
    case {'optical','optocondition'}
        TDT.SetTargetVal('ustim.StimDevice',1);
    case {'optoelectric','electrocondition'}
        TDT.SetTargetVal('ustim.StimDevice',2);   % both of el & opt
    case {'conditioning'}
        if ismember(csnum,[7])
            TDT.SetTargetVal('ustim.StimDevice',0);
        elseif ismember(csnum,[8])
            TDT.SetTargetVal('ustim.StimDevice',1);
        elseif ismember(csnum,[9])
            TDT.SetTargetVal('ustim.StimDevice',2);   % both of el & opt
        else
            TDT.SetTargetVal('ustim.StimDevice',3);   % no output from el/opt devices
        end
end

% if get(handles.togglebutton_ampblank,'Value')   % If amplifier blank is set
if 0   % If amplifier blank is set
    TDT.SetTargetVal('ustim.BlankAmp',1);
    TDT.SetTargetVal('ustim.BlankExtra',100);
else
    TDT.SetTargetVal('ustim.BlankAmp',0);
end












