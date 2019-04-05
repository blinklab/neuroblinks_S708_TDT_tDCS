clear
cd D:\My_docs\MATLAB

% --- session name ---
tgt_spk_file='LA001_171028_s01';


% loading mat file from spike2
load D:\My_docs\MATLAB\LA001_171028_s01


% loading video foleder
preprocess_1_extract_online_eye_2


data.marks.TrlN.valid


tr=8;

figure
subplot(211)
plot(trials.eye(tr).time,trials.eye(tr).trace)
xlim([-200 1600])

subplot(212)
plot(data.spkwf.CylP_time, double(data.spkwf.CylP(:,tr)))
xlim([-200 1600])


