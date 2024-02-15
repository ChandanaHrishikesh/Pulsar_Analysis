function oomp = dB_oomp(number)
% order of magnitude of number as power of 10: dB_oomp(792864) = 5, dB_oomp(0.00194752) = -3
oomp = floor(log10(abs(number)));