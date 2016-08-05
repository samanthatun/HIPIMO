function [ C ] = concO2( x,y,t )
%C Summary of this function goes here
%   Detailed explanation goes here

mkdir('diffusionfig1');
fig = figure('Color',[0.5 0.5 0.5]);

xmin = -100;
xmax = -xmin;
ymin = xmin;
ymax = xmax;

h= 2;
Nx = 1 + round((xmax-xmin)/h);
Ny = 1 + round((ymax-ymin)/h);

dt = 0.3;   %time step
D= 5.7;       %diffusion coefficient



C = zeros(Nx,Ny);
% C(11,18) = 15;   %source of C
% C(18,13) = 15;


for iter =1:100
%     C(18,18) = 15;
%     C(11,13) = 15;

    C(:,1) = 1;
    C(:,100) = 1;
    for i=2:Nx-1
        for j=2:Ny-1
            C(i,j) = C(i,j) + ((D*dt)/(h*h))*...
                (C(i-1,j)+C(i+1,j)+C(i,j-1)+C(i,j+1)- 4*C(i,j));
        end
    end

end
contourf(xmin:h:xmax, ymin:h:ymax, C,[0:0.1:2], 'edgecolor', 'none');
    axis equal;
    colormap(hsv)
    colorbar;
    pause(0.01);
    %fileName=['diffusionfig1/fig_',num2str(iter)];
    %saveas (fig, fileName, 'jpg');
end %end of loop for 1:100
