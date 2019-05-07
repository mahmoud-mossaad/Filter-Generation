function [ Filtered_signal ] = Digital_filter( zeros , poles ) % coord. of z/p chosen

Original_signal = load('ECG.mat');
%z = [1 2 -1 ]; % chosen zeroes
p1=poly(zeros); 
%p = 1; %chosen poles
p2=poly(poles);
Filtered_signal = filter(p1,p2,Original_signal.val); % p1 and p2 are coefficients of transfer function 
plot(Original_signal.val);
hold on
plot(Filtered_signal);



end

