data_folder = 'Z:\Carlos\ephys\video';
day = datestr(now,'yymmdd');
mice = {'CS004'};

for mouse=mice
    currentFolder = fullfile(data_folder,mouse{:},day);
    makeCompressedVideos(currentFolder,1);
    processSessions(currentFolder);
end
