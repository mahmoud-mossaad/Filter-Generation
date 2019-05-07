function [ gain ] = Get_gain_manual( zeros, poles, unit_circle_vector )
%unit_vector = [1; 1+i; 1;2+i;0;i]; 
%zeros =[1;i;0;2+i];
%poles = [ 1;2;i;2+2i];
den=1;
num = 1;

for i = 1:size(zeros)
    Distance(i) = norm (unit_circle_vector - zeros(i));
    num = Distance(i)* num;
    
end 
display(num)

for k = 1:size(poles)
    Distance(k) = norm (unit_circle_vector - poles(k));
    den = Distance(k)* den;
    
end 
display(den)

gain = num/den


end

