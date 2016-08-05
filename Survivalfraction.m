function [ S ] = Survivalfraction(x,y,t)
%SURVIVALFRACTION Summary of this function goes here
%   Detailed explanation goes here


d = 2;
n = 1;


    S = exp(-n*((Alpha(x,y,t)*d)+(Beta(x,y,t)*d*d)));

end

