function [ Y2 ] = betaOER( x, y, t )
%BETAOER Summary of this function goes here
%   Detailed explanation goes here

m1 = 3.25;
m2 = 1;
K = 3.28;
x = OxygenDynamics(x,y,t);


Y2 = ((K*(m1 - m2))/(x+K)) + m2;


end



