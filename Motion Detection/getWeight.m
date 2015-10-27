function w = getWeight( rows, cols )
%GETDISTANCE Summary of this function goes here
%   Detailed explanation goes here

center = [(rows-1)/2+1,(cols-1)/2+1];
d = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        d(j,i) = sqrt((i - center(2))^2 + (j - center(1))^2);
    end
end
d = d/(rows*cols);
w = 1 - d.^2;

end

