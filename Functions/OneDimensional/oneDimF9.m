function [ outVal, outGr ] = oneDimF9( x )
%Function f9 from functions.txt

outVal = -(x - 1)*(x - 8)*(x - 7)*(x - 8)*(x - 6);
outGr = -4192 + 3636*x - 1029*x^2 + 120*x^3 - 5*x^4;

end

