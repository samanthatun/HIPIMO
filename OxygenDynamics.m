function [ O2 ] = OxygenDynamics( x, y, t )
%Models oxygen in the tumor's environment
%   Detailed explanation goes here

%     if t <= 100
%         O2 = 1;
%     else O2 = .5;
        
    b = 10;
    a = .5;
    O2 = -a*x+b;


end

