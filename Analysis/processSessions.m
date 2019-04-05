function processSessions(folder)
% Traverse directories starting at ROOTFOLDER and recursively check for compressed video subdirectory, running processTrials when found.

TEST = 0;


% Check current directory for subdirectory containing the compressed files and if there is none, recursively go to the next nested directory
dirInfo= dir(folder);

isDir=[dirInfo.isdir];
dirNames={dirInfo(isDir).name};
dirNames(strcmp(dirNames, '.') | strcmp(dirNames, '..'))=[];

if isempty(dirNames)
	return	% This will allow us to stop recursion
end

for i=1:length(dirNames)
	if strcmp(dirNames{i},'compressed')
		disp(['Found folder ' fullfile(folder,dirNames{i})])
		if ~TEST && ~exist(fullfile(folder, 'trialdata.mat'),'file')
			% this is where we do what it is we want to do
			trials = processTrials(fullfile(folder,dirNames{i}));
			save(fullfile(folder, 'trialdata.mat'),'trials')
		end
	else
		processSessions(fullfile(folder,dirNames{i}))
	end

end

