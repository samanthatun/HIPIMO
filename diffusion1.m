
clear all;
clc;
close all;

mkdir('diffusionfig1');
fig = figure('Color',[0.5 0.5 0.5]);

xmin = -30;
xmax = -xmin;
ymin = xmin;
ymax = xmax;

h= 2;
Nx = 1 + round((xmax-xmin)/h);
Ny = 1 + round((ymax-ymin)/h);

dt = 0.3;   %time step
D= 2;       %diffusion coefficient

O2 = zeros(Nx,Ny);
% O2(11,18) = 15;   %source of O2
% O2(18,13) = 15;

for iter =1:250
%     O2(18,18) = 15;
%     O2(11,13) = 15;

    O2(:,1) = 5;


    for i=2:Nx-1
        for j=2:Ny-1
            O2(i,j) = O2(i,j) + ((D*dt)/(h*h))*...
                (O2(i-1,j)+O2(i+1,j)+O2(i,j-1)+O2(i,j+1)- 4*O2(i,j));
        end
    end
    
    contourf(xmin:h:xmax, ymin:h:ymax, O2,[0:0.1:2], 'edgecolor', 'none');
    axis equal;
    colormap(hsv)
    colorbar;
    pause(0.01);
    %fileName=['diffusionfig1/fig_',num2str(iter)];
    %saveas (fig, fileName, 'jpg');
end %end of loop for 1:100