function trials=processTrials(folder)
% TRIALS=processConditioningTrials(FOLDER)
% Return trials structure containing eyelid data and trial parameters for all trials in a session
% Optionally, specify threshold for binary image and maximum number of video frames per trial to use for extracting eyelid trace
% Note that this is a simplification of prior versions that required calibration info as an argument so it's not compatible with 
% function calls that assume additional arguments. This version should be more universally useful because it uses the calibration
% info stored in the metadata.


% if length(varargin) > 0
% 	thresh=varargin{1};
% end

% % This part is sloppy but I don't want to change all the other logic right now
% USESAVEDCALIB = 0;

% % Error checking
% if isstruct(calib)
% 	if ~isfield(calib,'scale') || ~isfield(calib,'offset')
% 		error('You must specify a valid calibration structure or file from which the structure can be computed')
% 	end
% elseif exist(calib,'file')
% 	[data,metadata]=loadCompressed(calib);
% 	if ~exist('thresh')
% 		thresh=metadata.cam.thresh;
% 	end
% 	[y,t]=vid2eyetrace(data,metadata,thresh,5);
% 	calib=getcalib(y);
% elseif isempty(calib)
% 	% use the calibration stored in each file's metadata
% 	USESAVEDCALIB = 1;
% else
% 	error('You must specify a valid calibration structure or file from which the structure can be computed')
% end

% By now we should have a valid calib structure to use for calibrating all files

if ~exist(folder,'dir')
	error('The directory you specified (%s) does not exist',folder);
end

% Get our directory listing, assuming the only MP4 files containined in the directory are the trials
% Later we will sort out those that aren't type='conditioning' based on metadata
fnames=getFullFileNames(folder,dir(fullfile(folder,'*.avi')));
% fnames=getFullFileNames(folder,dir(fullfile(folder,'*.mp4')));

% Only keep the files that match the pattern MOUSEID_DATE_SXX or MOUSEID_DATE_TXX, skipping for instance trials from Camera 2
% matches = regexp(fnames,'[A-Z]\d\d\d_\d\d\d\d\d\d_[ts]\d\d[a-z]?_\d\d\d','start','once');
% fnames = fnames(cellfun(@(a) ~isempty(a), matches));

% Preallocate variables so we can use parfor loop to process the files
eyelidpos=cell(length(fnames),1);	% We have to use a cell array because trials may have different lengths
tm=cell(length(fnames),1);			% Same for time

c_isi=NaN(length(fnames),1);
c_csnum=NaN(length(fnames),1);
c_csdur=NaN(length(fnames),1);
c_usdur=NaN(length(fnames),1);

e_depth=NaN(length(fnames),1);
e_freq=NaN(length(fnames),1);
e_pulsewidth=NaN(length(fnames),1);
e_traindur=NaN(length(fnames),1);
e_amp=NaN(length(fnames),1);
e_delay=NaN(length(fnames),1);

l_depth=NaN(length(fnames),1);
l_freq=NaN(length(fnames),1);
l_pulsewidth=NaN(length(fnames),1);
l_traindur=NaN(length(fnames),1);
l_amp=NaN(length(fnames),1);
l_delay=NaN(length(fnames),1);

t_ONOFF=NaN(length(fnames),1);
t_polarity=NaN(length(fnames),1);
t_amp=NaN(length(fnames),1);

trialnum=zeros(length(fnames),1);
ttype=cell(length(fnames),1);

trialtime=zeros(length(fnames),1);

numframes=zeros(length(fnames),1);

% Use a parallel for loop to speed things up
% matlabpool open	% Start a parallel computing pool using default number of labs (usually 4-8).
% cleaner = onCleanup(@() matlabpool('close'));

parfor i=1:length(fnames)

	[p,basename,ext]=fileparts(fnames{i});

    try
        [data,metadata]=loadCompressed(fnames{i});
    catch
        disp(sprintf('Problem with file %s', fnames{i}))
    end

    calib = struct;
 
	calib.scale = metadata.cam.calib_scale;
	calib.offset = metadata.cam.calib_offset;

	thresh = metadata.cam.thresh; 
    
    [eyelidpos{i},tm{i}]=vid2eyetrace(data,metadata,thresh,5,calib);

	c_isi(i)=metadata.stim.c.isi;
	c_csnum(i)=metadata.stim.c.csnum;
	c_csdur(i)=metadata.stim.c.csdur;
	c_usdur(i)=metadata.stim.c.usdur;

	e_depth(i)=metadata.stim.e.depth;
	e_freq(i)=metadata.stim.e.freq;
	e_pulsewidth(i)=metadata.stim.e.pulsewidth;
	e_traindur(i)=metadata.stim.e.traindur;
	e_amp(i)=metadata.stim.e.amp;
	e_delay(i)=metadata.stim.e.delay;

	l_depth(i)=metadata.stim.l.depth;
	l_freq(i)=metadata.stim.l.freq;
	l_pulsewidth(i)=metadata.stim.l.pulsewidth;
	l_traindur(i)=metadata.stim.l.traindur;
	l_amp(i)=metadata.stim.l.amp;
	l_delay(i)=metadata.stim.l.delay;
    
%     if isfield(metadata.stim,'t')
        t_ONOFF(i) = metadata.stim.t.ONOFF;
        t_polarity(i) = metadata.stim.t.polarity;
        t_amp(i) = metadata.stim.t.amp;
%     end
	
	trialnum(i)=metadata.cam.trialnum;
	ttype{i}=metadata.stim.type;

	trialtime(i)=metadata.ts(2);

	numframes(i)=length(eyelidpos{i});

	fprintf('Processed file %s\n',basename)
	
end

disp('Done reading data')

% matlabpool close


MAXFRAMES=max(numframes);


% Now that we know how long each trial is turn the cell arrays into matrices

traces=NaN(length(fnames),MAXFRAMES);
times=NaN(length(fnames),MAXFRAMES);
try
	for i=1:length(fnames) 
		trace=eyelidpos{i}; 
		t=tm{i}; 
		en=length(trace); 
		if en > MAXFRAMES
			en=MAXFRAMES; 
		end 
		traces(i,1:en)=trace(1:en); 
		times(i,1:en)=t(1:en); 
	end
catch
    disp(i)
end

track_cell = regexp(fnames,'_t?(\d\d)[a-d]?_','tokens','once');
track_number = cellfun(@safestr2double, track_cell);

trials.tracknum = track_number;

session_cell = regexp(fnames,'_s(\d\d)[a-d]?_','tokens','once');
session_number = cellfun(@safestr2double, session_cell);

trials.session_of_day = session_number;

trials.eyelidpos=traces;
trials.tm=times;
trials.fnames=fnames;

trials.c_isi=c_isi;
trials.c_csnum=c_csnum;
trials.c_csdur=c_csdur;
trials.c_usdur=c_usdur;

trials.e_depth=e_depth;
trials.e_freq=e_freq;
trials.e_pulsewidth=e_pulsewidth;
trials.e_traindur=e_traindur;
trials.e_amp=e_amp;
trials.e_delay=e_delay;

trials.l_depth=l_depth;
trials.l_freq=l_freq;
trials.l_pulsewidth=l_pulsewidth;
trials.l_traindur=l_traindur;
trials.l_amp=l_amp;
trials.l_delay=l_delay;

trials.t_ONOFF = t_ONOFF;
trials.t_polarity = t_polarity;
trials.t_amp = t_amp;

trials.trialnum=trialnum;
trials.type=ttype;

trials.iti=diff([0; trialtime]);
trials.iti(trials.trialnum==1) = 0;

function db=safestr2double(str)

if ~isempty(str)
    db=str2double(str);
else
    db=0;
end