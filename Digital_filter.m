function [ Filtered_signal ] = Digital_filter( zeros , poles,signal ) % coord. of z/p chosen


%z = [1 2 -1 ]; % chosen zeroes
p1 = poly(zeros);
p2 = poly(poles);
Filtered_signal = filter(p1,p2,signal); % p1 and p2 are coefficients of transfer function 

%plot(Original_signal.val);

%plot(Filtered_signal);



end

