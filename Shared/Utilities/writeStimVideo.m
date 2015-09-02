function writeStimVideo(data,metadata,fname,varargin)
    % WRITESTIMVIDEO(DATA,METADATA,FNAME,{FRAMERATE,FORMAT})
    % Write stim video to avi file FNAME including stim marker boxes. Optionally specify framerate, otherwise the video will be written using
    % the same frame rate as it was captured at. You can also optionally specify the video format to use (Default is Motion JPEG AVI). If your 
    % computer and Matlab version support it, MPEG-4 may give better compression ratios.

if isempty(data)
    % Die gracefully
    warning('There is no data to write.')
    return
end

if length(varargin) >= 1
    fr=varargin{1};
else
    fr=metadata.cam.fps;
end

if length(varargin) >= 2
    fmt=varargin{2};
else
    fmt='Motion JPEG AVI';
end


[m,n,p,t]=size(data);

% Make stim timing arrays


switch lower(metadata.stim.type)
    case 'none'
        % do nothing for now
        stim=zeros(m,n,3,t,'uint8');
    case 'puff'
        stfrm=round((metadata.cam.time(1)+metadata.stim.p.puffdelay)./1000.*metadata.cam.fps);
        enfrm=round((metadata.cam.time(1)+metadata.stim.p.puffdur)./1000.*metadata.cam.fps);
        stim=zeros(m,n,3,t,'uint8');
        stim(400:449,450:499,2,stfrm:stfrm+enfrm)=255;
    case 'electrical'       
        stfrm=round((metadata.cam.time(1)+metadata.stim.e.delay)./1000.*metadata.cam.fps);
        enfrm=round(metadata.stim.e.traindur./1000.*metadata.cam.fps);
        stim=zeros(m,n,3,t,'uint8');
        stim(400:449,500:549,1,stfrm:stfrm+enfrm)=255;       
    case 'conditioning'
        stfrmc=round(metadata.cam.time(1)./1000.*metadata.cam.fps); % for CS
        enfrmc=round(metadata.stim.c.csdur./1000.*metadata.cam.fps);
        stfrmu=round((metadata.cam.time(1)+metadata.stim.c.isi)./1000.*metadata.cam.fps);
        enfrmu=round(metadata.stim.c.usdur/1000.*metadata.cam.fps); % for US
        stim=zeros(m,n,3,t,'uint8');
        stim(400:449,500:549,3,stfrmc:stfrmc+enfrmc)=255;
        stim(400:449,550:599,2,stfrmu:stfrmu+enfrmu)=255;
    case 'optical'      
        stfrm=round((metadata.cam.time(1)+metadata.stim.l.delay)./1000.*metadata.cam.fps);
        enfrm=round(metadata.stim.l.traindur./1000.*metadata.cam.fps);
        stim=zeros(m,n,3,t,'uint8');
        stim(400:449,550:599,3,stfrm:stfrm+enfrm)=255;
    case 'optoelectric'      
        stfrme=round((metadata.cam.time(1)+metadata.stim.e.delay)./1000.*metadata.cam.fps);
        enfrme=round(metadata.stim.e.traindur./1000.*metadata.cam.fps);      
        stfrml=round((metadata.cam.time(1)+metadata.stim.l.delay)./1000.*metadata.cam.fps);
        enfrml=round(metadata.stim.e.traindur./1000.*metadata.cam.fps);
        stim=zeros(m,n,3,t,'uint8');
        stim(400:449,500:549,1,stfrme:stfrme+enfrme)=255;
        stim(400:449,550:599,3,stfrml:stfrml+enfrml)=255;
    case 'optocondition'      
        stfrmc=round((metadata.cam.time(1)+metadata.stim.l.delay)./1000.*metadata.cam.fps);
        enfrmc=round(metadata.stim.l.traindur./1000.*metadata.cam.fps);
        stfrmu=round((metadata.cam.time(1)+metadata.stim.p.puffdelay)./1000.*metadata.cam.fps);
        enfrmu=round((metadata.cam.time(1)+metadata.stim.p.puffdur)./1000.*metadata.cam.fps);
        stim=zeros(m,n,3,t,'uint8');
        stim(400:449,500:549,3,stfrmc:stfrmc+enfrmc)=255;
        stim(400:449,550:599,2,stfrmu:stfrmu+enfrmu)=255;
end


for i=1:t 
    F(i).cdata=repmat(data(:,:,:,i),[1 1 3 1])+stim(:,:,:,i); 
    F(i).colormap=[]; 
end

writeObj=VideoWriter(fname,fmt);
set(writeObj,'FrameRate',fr);
open(writeObj);

writeVideo(writeObj,F);
close(writeObj);

% hm=figure;
% set(hm,'Name','Instant Replay - close window when done','Position',[150 150 n m]);
% movie(hm,F,1,20);