function [ Avg ] = O2Lattice( t )
%calculates the amount of O2 for every x value and y value
%   loop that incorporates the x and y of the lattice

    total = 0;
    for x = 1:200;
        for y = 1:200
            total = total + OxygenDynamics(x,y,t);
        end
    end
    
    Avg = total/40000;


end

