function [distance] = EuclideanDistance(pc, j, n)

[x1,y1]=ind2sub(n,j);
[x2,y2]=ind2sub(n,pc);
distance = sqrt((x1-x2)^2+(y1-y2)^2); 
