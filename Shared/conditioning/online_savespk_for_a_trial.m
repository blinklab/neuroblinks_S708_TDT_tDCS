function online_savespk_for_a_trial(tnum)

ghandles=getappdata(0,'ghandles'); 
handles = guidata(ghandles.onetrialanagui);

% -- getting the data name for spike ---
codes1=get(handles.listbox_snip,'String');
codevalues=get(handles.listbox_snip,'Value');
Event=codes1{codevalues}; % 'None', 'Spks', 'Spk2', or 'LFPs'

metadata=getappdata(0,'metadata');
trials=getappdata(0,'trials');
TTX=getappdata(0,'ttx');

% --- init ----
trials.spk(tnum).y=NaN;
trials.spk(tnum).ylim=[NaN NaN];
trials.spk(tnum).time=NaN;
trials.spk(tnum).ts_interval=NaN;

% If using multi channel spike data and there are no spikes on channel 1, this function will generate an error 
% Select 'None' in the listbox to avoid the error
if strcmpi(Event,'None')
    setappdata(0,'trials',trials);
    return,    
end

% ---- get data from TDT ----    
dt=1;
Channel=1;
tlim=[-100 400]+200; % time from TrlN

[trln_ts,trln_values]=TDTgetEventData(TTX,'TrlN',0,0,'ALL');

ind1=find(trln_values==tnum);
if ~isempty(ind1)
%     trln_ts, trln_values
    StartTime=trln_ts(ind1(end))+tlim(1)/1000; EndTime=trln_ts(ind1(end))+(tlim(2)+100)/1000;
    [y t ts_interval] = getdataATR (TTX, Event, Channel, StartTime, EndTime);
    tind_1=find(t<(trln_ts(ind1(end))+tlim(2)/1000));  
    if ~isempty(tind_1), 
        tind_2=1:dt:tind_1(end);
        y=y(tind_2)*1e6;
        ylim=[min(y) max(y)];
        % time [ms], y [uV]
        if ~sum(abs(ylim))==0
            trials.spk(tnum).y=y;
            trials.spk(tnum).ylim=ylim;
            trials.spk(tnum).time=(t(tind_2)-trln_ts(ind1(end)))*1000-metadata.cam.time(1);
            trials.spk(tnum).ts_interval=ts_interval*1000;
        end
    end
end

% --- save results to memory ----
setappdata(0,'trials',trials);



function [y t ts_interval] = getdataATR (TTX, Event, Channel, StartTime, EndTime)
%Get the data and the timestamps; data will have the data organized as a
%series of colums each representing a block, timestamps will be a row where
%each value represents the time stamp of the first value in each block
%Also, since ReadEventsV tends to round up, add a tenth of a second to beginning 
%of the interval if the start isn't zero

TR = TTX.GetValidTimeRangesV();
if EndTime > TR(2), EndTime = TR(2); end
N = TTX.ReadEventsV(1000000, Event, Channel, 0, StartTime - 0.2, EndTime, 'ALL');
if N == 0
    y = NaN; t = NaN; ts_interval = NaN;
    return
end
y = TTX.ParseEvV(0,N);
timestamps = TTX.ParseEvInfoV(0,N,6);

%Organize the data into a meaningful waveform, where each value has a timestamp
%First construct an array with a time for every sample
ts_interval = 1/TTX.ParseEvInfoV(0,1,9);

t = timestamps(1) + (0:numel(y)-1) .* ts_interval;

%next organize all the data into one row containing all samples
y = reshape(y, 1, numel(y));

%Now we'll trim the excess samples that lie outside the specified range
k = 1;
if StartTime ~= 0
    while (t(k) < StartTime)
        k = k+1;
    end
end

j = length(t);
if EndTime ~= 0
    while (t(j) > EndTime)
        j = j-1;
    end
end

y = y(k:j);
t = t(k:j);



