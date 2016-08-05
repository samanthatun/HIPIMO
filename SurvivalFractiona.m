function [ F ] = SurvivalFractiona( x,y,t )
%SURVIVALFRACTIONA Summary of this function goes here
%   Detailed explanation goes here

d = 2;
n = 1;


    F = exp(-n*((StemCella(x,y,t)*d)+(StemCellb(x,y,t)*d*d)));
end

