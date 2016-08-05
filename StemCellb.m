function [ B ] = StemCellb( x,y,t )
%STEMCELLB Summary of this function goes here
%   Detailed explanation goes here
b = 0.05;

B = b/StemCellbOER(x,y,t);

end

