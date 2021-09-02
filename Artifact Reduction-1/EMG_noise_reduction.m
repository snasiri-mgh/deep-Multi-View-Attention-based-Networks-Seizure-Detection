%% This function takes 3-sec of 22 channel EEG signal, sampling frequency, and channels' labels as input. The outputs are reconstructed signal after movement, and EMG artifact reduction
% and deleted components.
function [rec,del_comp] = EMG_noise_reduction(eeg_sig,sfreq,hdrlabel) 
    eeg_sig1 = montaged(eeg_sig,hdrlabel);
    Xc = centre(eeg_sig1) ; % Centering signal
    [Scca,Acca]= EMGsorsep(Xc) ;
    Scca = pinv(Acca)*eeg_sig1 ; %computing CCA sources.
    
    A_max = max(abs(Acca)); 
    Scca_new = Scca.*A_max' ; % sources of CCA in signal domain.
    comp_to_del = find(max(abs(Scca_new')) > 200); % finding components related to movement artifact
    % EMG artifact reduction
    for i = 1 : 22
        s1(i) = Xcorr_emg(Scca_new(i,:),sfreq);
        s2(i) = bandpower(Scca_new(i,:),sfreq,[30,80])/ bandpower(Scca_new(i,:),sfreq,[1,80]);
    end
mm1 = find(s1 < 0.3 );
mm2 = find(s2 > 0.4 ) ;
q1 = zeros(22,1);
q2 = zeros(22,1);
q1(mm1)=1;
q2(mm2)=1;
sus_comp = find(q1.*q2 ==1);
for i = 1 : length(sus_comp) 
    [imf,residual] = emd(Scca(sus_comp(i),:));
    imf(:,[1,2])=0;
    Scca(sus_comp(i),:) = (sum(imf')' + residual)' ;
end
for i = 1 : 22
      s1(i) = Xcorr_emg(Scca(i,:),sfreq);
      s2(i) = bandpower(Scca(i,:),sfreq,[30,80])/ bandpower(Scca(i,:),sfreq,[1,80]);
end
mm1 = find(s1 < 0.3 );
mm2 = find(s2 > 0.4 ) ;
q1 = zeros(22,1);
q2 = zeros(22,1);
q1(mm1)=1;
q2(mm2)=1;
del_comp = find(q1.*q2 ==1);
Scca(del_comp,:) = 0;
Scca(comp_to_del,:) = 0;
rec =Acca*Scca ;
end

