function [mean_slope] = Xcorr_emg(input,sfreq)
x_corr = xcorr(input);
slope11 = abs(x_corr(3*sfreq + 2)/x_corr(3*sfreq +1));
slope12 = abs(x_corr(3*sfreq + 3)/x_corr(3*sfreq +1));
slope13 = abs(x_corr(3*sfreq + 4)/x_corr(3*sfreq +1));
slope14 = abs(x_corr(3*sfreq + 5)/x_corr(3*sfreq +1));

slope21 = abs(x_corr(3*sfreq - 2)/x_corr(3*sfreq +1));
slope22 = abs(x_corr(3*sfreq - 3)/x_corr(3*sfreq +1));
slope23 = abs(x_corr(3*sfreq - 4)/x_corr(3*sfreq +1));
slope24 = abs(x_corr(3*sfreq - 5)/x_corr(3*sfreq +1));

mean_slope = mean([slope11,slope12,slope13,slope14,slope21,slope22,slope23,slope24]);


end

