sfreq = hdr.samples(1);      %sampling frequency
%% EMG
[rec_emg,del_comp] = EMG_noise_reduction(eeg_sig,sfreq,hdr.label) ;    %EMG artifact reduction
% it takes eeg_sig ( 3-sec of recorded eeg), sampling frequency, and name
% of channels. rec_emg is reconstructed signal after emg artifact
% reduction.
%% EKG
[rec_emg_ekg] = ekg_artifact(ekg_channel,rec_emg,sfreq);  % This step takes reconstructed signal from
% previous step and remove EKG artifact. First input is ECG recording of
% signal. second one is reconstructed signal from previous step, and the
% third one is sampling frequency.
%% Eye
[reconstructed_signal,deleted_comp] = removing_eye_blink_artifact(rec_emg_ekg,sfreq);
% This step has an input of previous step and sampling frequency.