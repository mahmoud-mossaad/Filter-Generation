function [ gain ] = Get_gain_manual( zeros, poles, unit_circle_vector )
%unit_vector = [1; 1+i; 1;2+i;0;i]; 
%zeros =[1;i;0;2+i];
%poles = [ 1;2;i;2+2i];
den=1;
num = 1;
Distance = ones(length(unit_circle_vector),1);

for i = 1:size(unit_circle_vector)
    for j=1:size(zeros)
    Distance(i) = Distance(i)* norm (unit_circle_vector(i) - zeros(j));
    end
end 
display(num)

for i = 1:size(unit_circle_vector)
    for j=1:size(poles)
    Distance(i) = Distance(i)* 1/norm (unit_circle_vector(i) - poles(j));
    end
end 
plot(20*log(Distance));


end

