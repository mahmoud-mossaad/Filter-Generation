function [ gain ] = Get_gain_manual( zeros, poles,unit_circle_vector)

unit_circle_vector = load('data.mat');
unit_circle_vector = unit_circle_vector.data;
Distance = ones(length(unit_circle_vector),1);

for i = 1:size(unit_circle_vector)
    for j=1:size(zeros)
    Distance(i) = Distance(i)* norm (unit_circle_vector(i) - zeros(j));
    end
end 


for i = 1:size(unit_circle_vector)
    for j=1:size(poles)
    Distance(i) = Distance(i)* 1/norm (unit_circle_vector(i) - poles(j));
    end
end 
gain = 20*log10(Distance);


end

