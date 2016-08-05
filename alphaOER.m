function [ Y1 ] = alphaOER( x,y,t )
%OER Summary of this function goes here
%   Detailed explanation goes here

m1 = 1.75;
m2 = 1;
K = 3.28;





    Y1 = ((K*(m1 - m2))/(OxygenDynamics(x,y,t)+K)) + m2;



