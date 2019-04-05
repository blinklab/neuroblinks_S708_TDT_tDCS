%%% Goes through all folders/files in specified data directory & compresses
%%% the .mats. Assumes that already compressed files have been deleted,
%%% leaving only 4 figure files and the trialdata.mat file in an animal's
%%% folder.

parentdir='D:\Carlos';
cd(parentdir)
folders=dir('CS*');

for i=1:length(folders)

mouse = folders(i,1).name;
curdir=strcat('D:\Carlos\',mouse);
cd(curdir)
days=dir('17*');

for j=1:length(days)
    day=days(j,1).name;
    daydir=strcat(curdir,'\',day);
    cd(daydir)
    if exist('trialdata.mat')==0
        try
            trials = processTrials(fullfile(daydir,'compressed'),...
                fullfile(daydir,'compressed',sprintf('%s_%s_s01_calib.mp4',mouse,day)));
            save(fullfile(daydir, 'trialdata.mat'),'trials');
            [hf1,hf2, hf3, hf4]=makePlots(trials, 200, 1, 3, 7);
            
            hgsave(hf1,fullfile(cd,'CRs.fig'));
            hgsave(hf2,fullfile(cd,'CR_amp_trend.fig'));
            
            print(hf1,fullfile(cd,sprintf('%s_%s_CRs.pdf','CS002','170719')),'-dpdf')
            print(hf2,fullfile(cd,sprintf('%s_%s_CR_amp_trend.pdf','CS002','170719')),'-dpdf')
            close all
        catch ME
        end
    end

end

end