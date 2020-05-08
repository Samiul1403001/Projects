clear all;
clc;
A = xlsread('D:\Research_&_Projects\Blood pressure\data\WristCuff1.txt');
a1 = A(:,3);
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
windowSize = 5;
b = (1/windowSize)*ones(1,windowSize);
a=1;
f_a1 = filter(b,a,a1);
f_a1 = f_a1(20:length(a1));
figure;
plot(f_a1);

% Filtering
f_a2 = bandpass(f_a1, [.6 6.4], 20);
figure;
plot(f_a2);

% Slope determination
f_a3 = 100*diff(f_a2);
figure;
plot(f_a3);

% noise removal filter(low pass)
windowSize = 5;
b = (1/windowSize)*ones(1,windowSize);
a=1;
f_a4 = filter(b,a,f_a3);
f_a4 = f_a4(20:length(f_a4));
figure;
plot(f_a4);

% BP determination

