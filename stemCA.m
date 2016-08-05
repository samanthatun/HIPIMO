%%Base CA for stem - non-stem comparison chapter
%%30 DEC 2014

clear all;clf;close all;

figure('Position', [0 1600 600 500]) % make figure for CA display

% timestep counter
number = uicontrol('style','text', ...
    'string','1', ...
    'fontsize',12, ...
    'position',[20,400,50,20]);


%% simulation parameters
n  = 200; %domain size n x n
t_end = 1000; %time of simulation
age_limit = 5; % age limit of TACs 

%% cell behaviour
pSTEM_SymmDiv = .5; %probability of SYMMETRIC division of STEM CELLS
pSTEM_AsymmDiv = 0.5; %probability of ASYMMETRIC division STEM to TAC transition

pTAC_SymmDiv = 1; %probability of SYMMETRIC division of TAC CELLS
pTAC_Dediff = 0.001; %probability of TAC to STEM transition
pTAC_move = 0.2; %probability of motility of TAC CELLS

death_prob = 0.001; %frequency of random death
%% oxygen dynamics
    t_step = 10;
    tO = [0:t_step:t_end];
    p = zeros(length(tO),1);
    
    for i = 1:length(tO);
        p(i) = O2Lattice(tO(i));
    end
    
 %% linear oxygen progression
    
    t_step = .5;
    tG = [0:t_step:t_end];
    k = zeros(length(tG),1);
    
    for j = 1:length(tG);
        k(j) = OxygenDynamics(tG(j));
    end
    
%% OER

%     t_step = 1;
%     tK = [0:t_step:t_end];
%     r = zeros(length(tK),1);
%     
%     for d = 1:length(tK);
%         r(d) = OER(tK(d));
%     end
%     
%% Oxygen diffusion    
% mkdir('diffusionfig1');
% fig = figure('Color',[0.5 0.5 0.5]);
% 
% xmin = -30;
% xmax = -xmin;
% ymin = xmin;
% ymax = xmax;
% 
% h= 2;
% Nx = 1 + round((xmax-xmin)/h);
% Ny = 1 + round((ymax-ymin)/h);
% 
% dt = 0.3;   %time step
% D= 2;       %diffusion coefficient
% 
% drug = zeros(Nx,Ny);
% % drug(11,18) = 15;   %source of drug
% % drug(18,13) = 15;
% 
% for iter =1:250
% %     drug(18,18) = 15;
% %     drug(11,13) = 15;
% 
%     drug(:,1) = 5;
% 
% 
%     for i=2:Nx-1
%         for j=2:Ny-1
%             drug(i,j) = drug(i,j) + ((D*dt)/(h*h))*...
%                 (drug(i-1,j)+drug(i+1,j)+drug(i,j-1)+drug(i,j+1)- 4*drug(i,j));
%         end
%     end
% end
%     hold on
%     contourf(xmin:h:xmax, ymin:h:ymax, drug,[0:0.1:2], 'edgecolor', 'none');
%     axis equal;
%     colormap(hsv)
%     colorbar;
%     pause(0.01);
%     hold off
%% initialize the domain with all type 0 cells

cells = zeros(n,n);
TACage = cells; % initialize TACage with zeros

%% initialize domain with cells
cells(round(n/2),round(1+n/2)) = 0; %TAC CELLS
cells(round(1+n/2),round(n/2)) = 0; %TAC CELLS

cells(round(n/2),round(n/2)) = 0; %STEM CELLS
cells(round(1+n/2),round(1+n/2)) = 0.5; %STEM CELLS

%% initialize counters
totalcells = zeros(t_end,1);
totalTACcells = zeros(t_end,1);
totalSTEMcells = zeros(t_end,1);
TACprop = zeros(t_end,1);
STEMprop = zeros(t_end,1);
t = zeros(t_end,1);

cellsnew = cells; % create new matrix

%% main time loop
for j = 1:t_end
    
    update = randperm(n*n); % create random order list of all lattice sites
    
    for i = 1:n*n; % cycle through all lattice sites
        
        pc = update(i); % pick a random lattice site
        
        if cells(pc) == 0 % if the site is empty, skip loop
            continue
            
            %what to do if you pick a STEM CELL
        elseif cells(pc) == 0.5
            
            %find(cells(pc+randperm([-n n -1 +1])),0)  Jan's speed up
            [emptyp1,emptym1,emptypn,emptymn,emptypnp1,emptymnp1,emptymnm1,emptypnm1] = SpaceCheck8(pc,cellsnew,n); %check for space
            
            %if there is empty space, try to divide
            if emptyp1 || emptym1 || emptypn || emptymn || emptypnp1 || emptymnp1 || emptymnm1 || emptypnm1 == 1
                
                %roll die
                rn = rand;
                
                if rn <= pSTEM_SymmDiv
                    
                    %if you divide symmetrically
                    %choose one of the empty sites to fill
                    
                    rn = rand;
                    [cellsnew,TACage] = PlaceSTEMDaughter8(cellsnew,pc,rn,n,TACage);
                    
                elseif rn <= pSTEM_SymmDiv + pSTEM_AsymmDiv %else, divide asymmetrically
                    
                    rn = rand;
                    [cellsnew,TACage] = PlaceTACDaughter8(cellsnew,pc,rn,n,TACage);
                    
                end
            end
            
            
            
        else %what to do if you pick a TAC CELL
            
            %check for space
            [emptyp1,emptym1,emptypn,emptymn,emptypnp1,emptymnp1,emptymnm1,emptypnm1] = SpaceCheck8(pc,cellsnew,n);
            
            %if there is empty space, try to do something
            if emptyp1 || emptym1 || emptypn || emptymn || emptypnp1 || emptymnp1 || emptymnm1 || emptypnm1 == 1
                rn = rand;
                
                if rn <= pTAC_SymmDiv %probability of TAC cells growing
                    
                    %if you divide choose one of the empty sites to fill
                    rn = rand;
                    [cellsnew,TACage] = PlaceTACDaughter8(cellsnew,pc,rn,n,TACage);
                    
                elseif rn3 <= pTAC_SymmDiv + pTAC_move % TAC cell moving
                    
                    rn = rand;
                    [cellsnew] = Move(cellsnew,pc,rn,n); %move cell
                    
                elseif rn3 <= pTAC_SymmDiv + pTAC_move + pTAC_Dediff %TAC cell dedifferentiating
                    cellsnew(pc) = 0.5; %change to STEM
                    
                else cellsnew(pc)=cells(pc);
                    
                end
            end
        end
    end
    
    %% random and age related death
    for i=1:n*n
        if cellsnew(update(i)) == 1 && TACage(update(i)) > age_limit % kill TACs with age > TAC age limit
            cellsnew(update(i)) = 0;
            
        elseif cellsnew(update(i)) == 1 && death_prob > 0    %randomly kill some TAC cells at frequency defined above
            kill = rand;
            if kill < death_prob
                cellsnew(update(i)) = 0;
            end
        end
    end
    
    cells=cellsnew;
    
    %% visualise each X timesteps
    if mod(j,1)==0 && j > 1
        hold on
        vis = image(50*cellsnew);
        %imh = image(cat(3,cellsnew,z,z));
        %set(vis, 'erasemode', 'none')
        colormap jet
        axis equal
        axis tight
        hold off
        
        pause(0.0000000000001)
    end
    
    totalTACcells(j)=length(find(cellsnew==1.0));
    totalSTEMcells(j)=length(find(cellsnew==0.5));
    totalcells(j)=totalTACcells(j)+totalSTEMcells(j);
    TACprop(j)=totalTACcells(j)/totalcells(j);
    STEMprop(j)=totalSTEMcells(j)/totalcells(j);
    t(j)=j;
    
    if totalcells(j)==0 %|| min(min(cellsnew)) == 0.5
        %t_end = j;
        break
    else
    end
    
    stepnumber = 1 + str2double(get(number,'string'));
    set(number,'string',num2str(stepnumber))

x = 0:1:200;
figure(1)
filename = 'stemCAmovie.gif';

frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if j == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end  
end
 
% 
%  for   k = 1:100;
% 	plot(fft(eye(k+16)))
% 	axis([0 200 0 200])
% 	M(k) = getframe();
%  end
%  


% figure()
% plot(tO,p)
% % figure() 

plot_start = 1;

figure('Position', [600 1600 1100 500])

subplot(2,2,1)
h=plot(t(plot_start:end),totalSTEMcells(plot_start:end),t(plot_start:end),totalTACcells(plot_start:end),t(plot_start:end),totalcells(plot_start:end),'.');
set(h,'linewidth',3);
set(h,'markersize',10);
set(gca,'fontsize',16);
title({'TAC-STEM CA'});
legend('STEM','TAC','total');
xlabel('timestep');
ylabel({'proportion';'of population'});

subplot(2,2,2)
o=plot(t(plot_start:end),STEMprop(plot_start:end),t(plot_start:end),TACprop(plot_start:end));
set(o,'linewidth',3);
set(o,'markersize',10);
set(gca,'fontsize',16);
title({'proportions'});
legend('STEM','TAC');

% subplot(2,2,3)
% plot(tO,p,'linewidth', 2)
% 
% subplot(2,2,4)
% plot(tG,k,'linewidth', 2)

