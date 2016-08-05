function [ Avg ] = SFLattice2( t )
%SFLATTICE2 Summary of this function goes here
%   Detailed explanation goes here
total = 0;
    for x = 1:200;
        for y = 1:200
            total = total + SurvivalFractiona(x,y,t);
       
        end
    end
    
    Avg = total/40000;

end

