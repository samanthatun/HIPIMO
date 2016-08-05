function [ B ] = Beta(x,y,t)
%BETA Summary of this function goes here
%   Detailed explanation goes here

b = 0.05;

B = b/betaOER(x,y,t);
end

