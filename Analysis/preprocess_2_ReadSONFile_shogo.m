% function matfname=preprocess_2_ReadSONFile_shogo(filename,varargin)
clear
filename='LA001_171028_s01.smr';
nargin=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MATFILENAME = READSONFILE(SONFILENAME,{START,END})
% Create MAT file data structure from TDT tank data for optoelectrophys.
% Must run WriteSONFile() and sort spikes before running this
% script. 
% Optionally supply start and end times to read in seconds as additional arguments.

% % ----- for units --------
if exist(['S:\DATA_2012_upenn\units_smr\' filename],'file'),
% if str2num(filename(6:11))<151231,
    path1='S:\DATA_2012_upenn\units_smr\';
    savepath='S:\DATA_2012_upenn\units_mat\';
elseif exist(['S:\DATA_2012_upenn\units_unsorted\' filename],'file'),
    path1='S:\DATA_2012_upenn\units_unsorted\';
    savepath='S:\DATA_2012_upenn\units_mat\';
elseif exist(['S:\DATA_2012_upenn\bhv_conditioning\' filename],'file'),
    path1='S:\DATA_2012_upenn\bhv_conditioning\';
    savepath='S:\DATA_2012_upenn\bhv_mat\';
elseif exist(['S:\DATA_2015_baylor\units_unsorted\' filename],'file'),
    path1='S:\DATA_2015_baylor\units_unsorted\';
    savepath='S:\DATA_2015_baylor\units_mat\';
elseif exist(['S:\DATA_2015_baylor\units_smr\' filename],'file'),
    path1='S:\DATA_2015_baylor\units_smr\';
    savepath='S:\DATA_2015_baylor\units_mat\';
elseif exist(['D:\My_docs\MATLAB\' filename],'file'),
    path1='D:\My_docs\MATLAB\';
    savepath='D:\My_docs\MATLAB\';
%     savepath='S:\DATA_2015_baylor\units_mat\';
end

addpath('S:\ShareFolder\MATLAB_work\mfile_common');   addpath_shogo;

sonfname=[path1 filename];

switch nargin
    case 1
        st=0;                  en=0;
    case 2
        st=varargin{1};        en=0;
    case 3
        st=varargin{1};        en=varargin{2};
end

% Read spktime data from SMR
fid=fopen(sonfname);
smrchanlist=SONChanList(fid);

% ------- extract whole spike wave  %% important waveform should be the first channle (ch 16)-----
for i=1:length(smrchanlist)
    if strcmpi(smrchanlist(i).title(1:4),'LFPs'),
        [LFP,h_LFP]=SONGetChannel(fid,smrchanlist(i).number); % Read data
    end
    
    if str2num(filename(6:11))<151231, % for Penn
        if strcmpi(smrchanlist(i).title(1:4),'Spks'),
            [waveform,h_waveform]=SONGetChannel(fid,smrchanlist(i).number); % Read data
        elseif length(smrchanlist(i).title)>5, 
            if strcmpi(smrchanlist(i).title(1:6),'Unit01'),
                [waveform,h_waveform]=SONGetChannel(fid,smrchanlist(i).number); % Read data
            end
        end
%         if strcmpi(smrchanlist(i).title(1:4),'Spk2'),
%             [waveform,h_waveform]=SONGetChannel(fid,smrchanlist(i).number); % Read data
%         end
    else % for Baylor
        if strcmpi(smrchanlist(i).title(1:4),'Spk2'),
            [waveform,h_waveform]=SONGetChannel(fid,smrchanlist(i).number); % Read data
        end
        if strcmpi(smrchanlist(i).title(1:4),'m_Sp'),
            [waveform,h_waveform]=SONGetChannel(fid,smrchanlist(i).number); % Read data
        end
        
        if strcmpi(smrchanlist(i).title(1:4),'CylP'),
            [CylP,h_CylP]=SONGetChannel(fid,smrchanlist(i).number); % Read data
        end
    end
end
if exist('waveform','var'),  waveform=[waveform;NaN*ones(300,1)];  end

% Sort through the channels in SMR file and create MAT vars from channels
for i=1:length(smrchanlist)
    switch smrchanlist(i).kind
        case {2,3} % Events (time stamp of spikes)
            if smrchanlist(i).number < 17 % Spike times
                fieldname=sprintf('unit%02d',smrchanlist(i).number);
                [spkdata,h_spkdata]=SONGetChannel(fid,smrchanlist(i).number);
%                 spkdata(1)=[]; spkdata(end)=[];
                data.neuron.(fieldname).spiketimes=spkdata;
                
                %%%%%%% added by Shogo %%%%%%%%
                sampleinterval=h_waveform.sampleinterval;
                samplerate=1000000/sampleinterval;
                spktimings=spkdata;
                % --- detecting t=zero ---
                wf_win=[1:20]';
                spkadc=waveform((wf_win)*ones(1,length(spktimings))...
                    +ones(length(wf_win),1)*round(spktimings'*samplerate));
                
                med_wf=median(spkadc,2);
                [y,ind0]=max(double(abs(med_wf)).*(0.90.^[1:length(med_wf)]'));
                % --- end of detection for t=zero ---
                if smrchanlist(i).number<=10, % -- for SS
                    wf_win=[-6:24]';   remove_thr=1.5;
                else, % -- for CS
                    wf_win=[-128:256]'; remove_thr=0.15;
                end
                data.neuron.(fieldname).wave=waveform((wf_win+ind0)*ones(1,length(spktimings))...
                    +ones(length(wf_win),1)*round(spktimings'*samplerate));
                
                %                     data.neuron.(fieldname).wave=spkdata.adc; 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
%                 [events,h_events]=SONGetChannel(fid,smrchanlist(i).number);           
%                 chantitle=smrchanlist(i).title;
%                 if isempty(events) || max(events)==0
%                     continue
%                 end
%                 data.events.(chantitle).times=events;
%                 data.events.(chantitle).valid=ones(size(events));
            end          
        case 6 
            % Wavemark channels treated differently than event channels containing spike times 
            % since we might want to import spike shapes too (not implemented yet).
            if smrchanlist(i).number < 17 % Spike times
                fieldname=sprintf('unit%02d',smrchanlist(i).number);
                [spkdata,h_spkdata]=SONGetChannel(fid,smrchanlist(i).number);
                if isempty(spkdata)
                    continue
                end
                data.neuron.(fieldname).spiketimes=spkdata.timings;
                data.neuron.(fieldname).sortcode=spkdata.markers;
                
%                 %%%%%%% added by Shogo %%%%%%%%
%                 samplerate=1000000/h_spkdata.sampleinterval;
%                 spktimings=spkdata.timings;
%                 % --- detecting t=zero ---
%                 spkadc=spkdata.adc;
%                 med_wf=median(spkadc,2);
%                 [y,ind0]=max(double(abs(med_wf)).*(0.90.^[1:length(med_wf)]'));
%                 % --- end of detection for t=zero ---
%                 if smrchanlist(i).number<=10, % -- for SS
%                     wf_win=[-6:24]';   remove_thr=1.5;
%                 else, % -- for CS
%                     wf_win=[-24:244]'; remove_thr=0.15;
%                 end
%                 data.neuron.(fieldname).wave=waveform((wf_win+ind0)*ones(1,length(spktimings))+ones(length(wf_win),1)*round(spktimings'*samplerate));
%                 %                     data.neuron.(fieldname).wave=spkdata.adc; 
%                 
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
        case 8 % Textmarks
            chantitle=smrchanlist(i).title;
            [marks,h_marks]=SONGetChannel(fid,smrchanlist(i).number); % Read data
            % Generate field name and assign all data for channel to
            % dynamically generated field name
            
            
            if isempty(marks) || ~isstruct(marks)
                continue
            end
            values=zeros(h_marks.npoints,1);
            valuechar=char(marks.text');
            [m,n]=size(valuechar);
            for j=1:m 
                values(j)=str2num(valuechar(j,:)); 
            end
            % Generate struct field: s.(name) dynamically generates fieldname
            data.marks.(chantitle).times=marks.timings;
            data.marks.(chantitle).values=values;
            data.marks.(chantitle).valid=ones(size(marks.timings));
            
            %%%%%%% added by Shogo %%%%%%%%
            if strcmpi(chantitle,'TrlN'),
                data.marks.(chantitle).valid=double(marks.markers(:,1));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    end
end
fclose(fid);

% Save some metadata so we know which version of SMR file we used
data.date=datestr(now);
data.smrfilename=sonfname;

%%%%%%%%%% compensation for some data %%%%%%%%%
% --- LED channel correction only for T022 ----
if sum(ismember(str2double(filename(2:4)),[022])), % swap 1 and 2
    data.marks.CsCh.values(data.marks.CsCh.values==1)=-1;
    data.marks.CsCh.values(data.marks.CsCh.values==2)=1;
    data.marks.CsCh.values(data.marks.CsCh.values==-1)=2;
end
% --- LED channel correction only for T022 ----
if sum(ismember(str2double(filename(6:11)),[140902 140917]))
    [data] = laser_data_compensation(data);
end
%%%%%%%%%% end of compensation %%%%%%%%%%%%



% --- extract spike waveform data ---
spkwf_int=[-520 1220];

ntr=length(data.marks.TrlN.times);
valid_tr=data.marks.TrlN.valid>0;
timings=data.marks.TrlN.times(valid_tr);
% ------- data for each trial ---------
if exist('waveform','var'),
    spkwf_ind=round(spkwf_int/h_waveform.sampleinterval*1000);
    wf_win=[spkwf_ind(1):spkwf_ind(2)]';
    spkwf.trace=NaN*ones(length(wf_win),ntr);
    % --- extract data ---
    spkwf.trace(:,valid_tr)=waveform((wf_win)*ones(1,length(timings))...
        +ones(length(wf_win),1)*round(timings'*1000000/h_waveform.sampleinterval));
    spkwf.time=wf_win*h_waveform.sampleinterval/1000;
end
% --- extract LFP ----
if exist('LFP','var') & length(LFP)>10
LFP_whole_time=[1:length(LFP)]*h_LFP.sampleinterval*10^-3; % ms
spkwf.LFP_time=spkwf_int(1):1:spkwf_int(2);
spkwf.LFP=NaN*ones(length(spkwf.LFP_time),ntr);
for tr=find(valid_tr')
    tm_algn=LFP_whole_time-data.marks.TrlN.times(tr)*1000;
    tind=tm_algn>=spkwf_int(1)-2 & tm_algn<=spkwf_int(2)+2;
    spkwf_LFP(:,tr)=interp1(tm_algn(tind),double(LFP(tind)),spkwf.LFP_time);
end
spkwf.LFP=int16(spkwf_LFP*2);
spkwf.LFP_scale=0.30517578125/2;
end

% --- extract cylinder data ----
if exist('CylP','var')
    CylP=diff(CylP);
    CylP(CylP<-1e4)=20;  CylP(CylP>1e4)=-20;  CylP=sign(CylP);
    CylP_whole_time=[1:length(CylP)]*h_CylP.sampleinterval*10^-3; % ms
    spkwf.CylP_time=spkwf_int(1):1:spkwf_int(2);
    spkwf.CylP=NaN*ones(length(spkwf.CylP_time),ntr);
    for tr=find(valid_tr')
        tm_algn=CylP_whole_time-data.marks.TrlN.times(tr)*1000;
        tind=tm_algn>=spkwf_int(1)-2 & tm_algn<=spkwf_int(2)+2;
        spkwf_CylP(:,tr)=interp1(tm_algn(tind),double(CylP(tind)),spkwf.CylP_time, 'nearest');
    end
    spkwf.CylP=int8(spkwf_CylP);
    spkwf.CylP_times=CylP_whole_time/1000;
    spkwf.CylP_whole=CylP;
    % spkwf.CylP_scale=0.30517578125/2;
end

% spkwf_ind=round(spkwf_int/h_LFP.sampleinterval*1000);
% wf_win=[spkwf_ind(1):spkwf_ind(2)]';
% spkwf.LFP=NaN*ones(length(wf_win),ntr);
% spkwf.LFP(:,valid_tr)=LFP((wf_win)*ones(1,length(timings))...
%           +ones(length(wf_win),1)*round(timings'*1000000/h_LFP.sampleinterval));
%       
% spkwf.LFP_time=wf_win*h_LFP.sampleinterval/1000;

data.spkwf=spkwf;

% Save the data structure
% Subsequent analysis routines will leave data untouched (except maybe
% events.x.valid) and will save analysis results into a 'results' struct.
[path,basename,ext]=fileparts(sonfname);
matfname=[basename '.mat'];
save([savepath matfname],'data')






