function [ O2diff ] = Diffusion( x,y,t )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mkdir('diffusionfig1');
fig = figure('Color',[0.5 0.5 0.5]);

xmin = 0;
xmax = 200;
ymin = xmin;
ymax = xmax;

h= 2;
Nx = 1 + round((xmax-xmin)/h);
Ny = 1 + round((ymax-ymin)/h);

dt = 0.3;   %time step
D= 5.7;       %diffusion coefficient



O2diff = zeros(Nx,Ny);
% O2diff(11,18) = 15;   %source of O2diff
% O2diff(18,13) = 15;


for iter =1:10000
%     O2diff(18,18) = 15;
%     O2diff(11,13) = 15;

    O2diff(:,1) = 10;


    for i=2:Nx-1
        for j=2:Ny-1
            O2diff(i,j) = O2diff(i,j) + ((D*dt)/(h*h))*...
                (O2diff(i-1,j)+O2diff(i+1,j)+O2diff(i,j-1)+O2diff(i,j+1)- 4*O2diff(i,j));
        end
    end

end
contourf(xmin:h:xmax, ymin:h:ymax, O2diff,[0:0.1:2], 'edgecolor', 'none');
    axis equal;
    colormap(jet)
    colorbar;
    pause(0.01);
    %fileName=['diffusionfig1/fig_',num2str(iter)];
    %saveas (fig, fileName, 'jpg');
end %end of loop for 1:100
