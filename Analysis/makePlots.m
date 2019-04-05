function varargout = makePlots(trials,varargin)

if length(varargin) > 0
    isi = varargin{1};
    session = varargin{2};
    us = varargin{3};
    cs = varargin{4};
else
    isi = 200;
    session = 1; % Session s01
    us = 3;
    cs = 7;
end

%% Eyelid traces
hf1=figure;
hax=axes;

%idx = find(trials.c_usnum==us & trials.c_csnum==cs & ismember(trials.session_of_day,session));
idx = find(trials.c_usnum==us & trials.c_csnum==cs & ismember(trials.session_of_day,session));

set(hax,'ColorOrder',jet(length(idx)),'NextPlot','ReplaceChildren');
plot(trials.tm(1,:),trials.eyelidpos(idx,:)')
hold on
plot(trials.tm(1,:),mean(trials.eyelidpos(idx,:)),'k','LineWidth',2)
axis([trials.tm(1,1) trials.tm(1,end) -0.1 1.1])
title('CS-US')
xlabel('Time from CS (s)')
ylabel('Eyelid pos (FEC)')


%% CR amplitudes
pre = 1:tm2frm(0.1);
win = tm2frm(0.2+isi/1e3):tm2frm(0.2+isi/1e3+0.015);
cramp = mean(trials.eyelidpos(:,win),2) - mean(trials.eyelidpos(:,pre),2);

hf2=figure;

plot(trials.trialnum(idx),cramp(idx),'.')
hold on
plot([1 length(trials.trialnum)],[0.1 0.1],':k')
axis([1 length(trials.trialnum) -0.1 1.1])
title('CS-US')
xlabel('Trials')
ylabel('CR size')

%% CR eyelid traces
hf3 = figure;
hax2 = axes;

set(hax2,'ColorOrder',jet(length(idx)),'NextPlot','ReplaceChildren');
plot(trials.tm(1,73:81), trials.eyelidpos(idx,73:81));
axis([trials.tm(1,73) trials.tm(1,81) -0.1 1.1]);
title('CR 40ms before US');
xlabel('Time from CS (s)');
ylabel('Eyelid pos (FEC)');


%% UR velocity
hf4 = figure;

v = diff(mean(trials.eyelidpos(idx, :)));
plot(trials.tm(1,1:199), diff(trials.eyelidpos(idx,:),1,2));
hold on;
plot(trials.tm(1,1:199), v, 'k', 'Linewidth', 2);

title('CR Velocity');
xlabel('Time from CS (s)');
ylabel('Velocity (FEC/5ms)');


%%
if nargout > 0
    varargout{1}=hf1;
    varargout{2}=hf2;
    varargout{3}=hf3;
    varargout{4}=hf4;
end