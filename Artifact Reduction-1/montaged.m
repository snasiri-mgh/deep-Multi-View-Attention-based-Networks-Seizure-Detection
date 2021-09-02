function [montage] = montaged(E,hdrlabel)
Fp1 = find(hdrlabel == ["EEGFP1REF"]);
Fp2 = find(hdrlabel == ["EEGFP2REF"]);
F3 = find(hdrlabel == ["EEGF3REF"]);
F4 = find(hdrlabel == ["EEGF4REF"]);
C3 = find(hdrlabel == ["EEGC3REF"]);
C4 = find(hdrlabel == ["EEGC4REF"]);
P3 = find(hdrlabel == ["EEGP3REF"]);
P4 = find(hdrlabel == ["EEGP4REF"]);
O1 = find(hdrlabel == ["EEGO1REF"]);
O2 = find(hdrlabel == ["EEGO2REF"]);
F7 = find(hdrlabel == ["EEGF7REF"]);
F8 = find(hdrlabel == ["EEGF8REF"]);
T3 = find(hdrlabel == ["EEGT3REF"]);
T4 = find(hdrlabel == ["EEGT4REF"]);
T5 = find(hdrlabel == ["EEGT5REF"]);
T6 = find(hdrlabel == ["EEGT6REF"]);
A1 = find(hdrlabel == ["EEGA1REF"]);
A2 = find(hdrlabel == ["EEGA2REF"]);
Fz = find(hdrlabel == ["EEGFZREF"]);
Cz = find(hdrlabel == ["EEGCZREF"]);
Pz = find(hdrlabel == ["EEGPZREF"]);
montage0 = E(Fp1,:) - E(F7,:);
montage1 = E(F7,:) - E(T3,:);
montage2 = E(T3,:) - E(T5,:);
montage3 = E(T5,:) - E(O1,:);
montage4 = E(Fp2,:) - E(F8,:);
montage5 = E(F8,:) - E(T4,:);
montage6 = E(T4,:) - E(T6,:);
montage7 = E(T6,:) - E(O2,:);
montage8 = E(A1,:) - E(T3,:);
montage9 = E(T3,:) - E(C3,:);
montage10 = E(C3,:) - E(Cz,:);
montage11 = E(Cz,:) - E(C4,:);
montage12 = E(C4,:) - E(T4,:);
montage13 = E(T4,:) - E(A2,:);
montage14 = E(Fp1,:) - E(F3,:);
montage15 = E(F3,:) - E(C3,:);
montage16 = E(C3,:) - E(P3,:);
montage17 = E(P3,:) - E(O1,:);
montage18 = E(Fp2,:) - E(F4,:);
montage19 = E(F4,:) - E(C4,:);
montage20 = E(C4,:) - E(P4,:);
montage21 = E(P4,:) - E(O2,:);






montage = [montage0;montage1;montage2;montage3;montage4;montage5;montage6;montage7;montage8;montage9...
    ;montage10;montage11;montage12;montage13;montage14;montage15;montage16;montage17;montage18;montage19...
    ;montage20;montage21];
end

