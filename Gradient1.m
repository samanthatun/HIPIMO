
clear all;
clc;
close all;

%mkdir('gradient');
fig = figure('Color',[0.5 0.5 0.5]);

xmin = 0;
xmax = 500;
ymin = xmin;
ymax = xmax;

h= 2;
Nx = 1 + round((xmax-xmin)/h);
Ny = 1 + round((ymax-ymin)/h);

dt = 1;   %time step
D= .5;       %diffusion coefficient

drug = zeros(Nx,Ny);

%plot (0, -20,'ko', 'MarkerSize', 50)
% x = [0;20;40;-20;-40]
% y = [0;0;0;0;0]
% scatter3(x,y,10*ones(size(x)), 'w')

% rad = .25;
% dt = .5;
% for i = 1:100
% 
%     numCells = size (x, 1)
% 
%     for j=1:numCells
%        ang = rand()*2*pi;      %cell moves
%        x(j) = x(j) + rad*[cos(ang)]*dt;
%        y(j) = y(j) + rad*[ sin(ang)]*dt;
%     end
%     pause(.001);
% end 

for iter =1:1000
        

    for k = 1:Nx
        drug(Ny,k)= 0; 
        drug(1,k) = 1;
        drug(k,1) = drug(k,2); 
        drug(k,Nx) = drug(k,Nx-1); 
    end 
    

    for i=2:Nx-1
        for j=2:Ny-1
            drug(i,j) = drug(i,j) + ((D*dt)/(h*h))*...
                (drug(i-1,j)+drug(i+1,j)+drug(i,j-1)+drug(i,j+1)- 4*drug(i,j));
        end
    end
    
    if mod(iter,5000) == 0
        %contourf(xmin:h:xmax, ymin:h:ymax, drug',[0:0.05:1], 'edgecolor', 'none');
        [A,B] = meshgrid(xmin:h:xmax, ymin:h:ymax);
        surface(A,B,drug','edgecolor', 'none');
        xlim([xmin xmax])
        ylim([ymin ymax])
        axis equal;
        colormap(jet)
        colorbar;
        %pause
        view([0 90])
        pause(0.05);
       
    end
    
    %fileName=['gradient/fig_',num2str(iter)];
    %saveas (fig, fileName, 'jpg');
       %fname = 'gradient.txt';
       
  
end   
size(drug)

save('gradient.txt','drug','-ascii');