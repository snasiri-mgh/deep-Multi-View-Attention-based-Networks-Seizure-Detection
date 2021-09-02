%% function of eye artifact reduction
function [reconstructed_signal,deleted_comp] = removing_eye_blink_artifact(input_signal_montage,sfreq)
    P=15;
    Xc = centre(input_signal_montage) ;
    [Acom2] = COM2R(Xc,P) ;
    Scom2 = pinv(Acom2)*input_signal_montage ;
    bp_eye = bandpower(Scom2(1:5,:)',sfreq,[1,2])./bandpower(Scom2(1:5,:)',sfreq,[1,30]);
    first_step_eye_source = find(bp_eye > 0.2);
    channel_amp = sum(abs(Acom2([1,5,15,19],1:5)))./sum(abs(Acom2(:,1:5)));
    channel_amp_selected_channel = find(channel_amp > 0.4);
    m = zeros(5,1);
    n = zeros(5,1);
    m(first_step_eye_source) = 1;
    n(channel_amp_selected_channel) = 1;
    total_selected_comp = find(m.*n == 1);
    deleted_comp = Scom2(total_selected_comp,:);
    Scom2(total_selected_comp,:) = 0;
    reconstructed_signal = Acom2*Scom2 ;
end

