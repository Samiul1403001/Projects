clear all;
clc;
A = xlsread('data.xlsx');
a1 = A(:,2);
a1(sum(isnan(a1), 2) == 1, :) = [];
[a1max,Pa1max] = max(a1);
a1 = a1(Pa1max:length(a1));
figure;
plot(a1);
a1_diff = diff(a1);
[a1_diff_max, Pa1_diff_max] = max(a1_diff);
a1 = a1(Pa1max:Pa1_diff_max-1);
figure;
plot(a1);

% noise removal filter(low pass)
windowSize = 10;
b = (1/windowSize)*ones(1,windowSize);
a=1;
f_a1 = filter(b,a,a1);
f_a1 = f_a1(20:length(a1));
figure;
plot(f_a1);

% Slope determination
f_a2 = diff(f_a1);
figure;
plot(f_a2);

% noise removal filter(low pass)
windowSize = 10;
b = (1/windowSize)*ones(1,windowSize);
a=1;
f_a3 = filter(b,a,f_a2);
figure;
plot(f_a3);

% Filtering
f_a4 = bandpass(f_a3, [.6 6.4], 20);
figure;
plot(f_a4);

% BP determination

