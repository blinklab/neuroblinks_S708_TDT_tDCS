function InitCam(ch)

% First delete any existing image acquisition objects

imaqreset

% vidobj = videoinput('avtmatlabadaptor64_r2009b', 1, 'Mono8_640x480_Binning_1x1');
vidobj = videoinput('gige', ch, 'Mono8');
src = getselectedsource(vidobj);
src.ExposureTimeAbs = 4900;
% src.AllGainRaw=7;
% src.AllGainRaw=12;
src.GainRaw=12;		% The name changed in Matlab 2013a
% src.NetworkPacketSize = '9014';
% src.PacketSize = '9014';
% src.PacketSize = 9014;  % Matlab 2013a uses integer value
src.PacketSize = 8228;	% Matlab 2013a uses integer value
src.StreamBytesPerSecond=124e6; % Set based on AVT's suggestion
vidobj.LoggingMode = 'memory'; vidobj.FramesPerTrigger=1;
src.AcquisitionFrameRateAbs=200;

%for 2010b
% triggerconfig(vidobj,'Hardware','RisingEdge','externalTrigger')
% triggerconfig(vidobj, 'hardware', 'risingEdge', 'externalTrigger');
% triggerconfig(vidobj, 'hardware', 'RisingEdge', 'Line1-FrameStart');

%for 2013a
triggerconfig(vidobj, 'hardware', 'DeviceSpecific', 'DeviceSpecific');
src.TriggerMode='On';
src.TriggerSelector='FrameStart';
src.TriggerSource='Freerun';
% src.TriggerSource='Line1';

%% Save objects to root app data

setappdata(0,'vidobj',vidobj)
setappdata(0,'src',src)