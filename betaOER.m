function [ Y2 ] = betaOER( x,y,t )
%BETAOER Summary of this function goes here
%   Detailed explanation goes here

m1 = 3.25;
m2 = 1;
K = 3.28;



Y2 = ((K*(m1 - m2))/((OxygenDynamics(x,y,t))^2+K)) + m2;


end



