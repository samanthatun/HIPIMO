function [cellsnew,TACage] = PlaceSTEMDaughter8(cellsnew ,x,rn,n,TACage)
% this function moves a cell to a random empty neighbouring location

[emptyp1,emptym1,emptypn,emptymn,emptypnp1,emptymnp1,emptymnm1,emptypnm1] = SpaceCheck8(x,cellsnew,n); %check for space
naybur = emptyp1 + emptym1 + emptypn + emptymn + emptypnp1 + emptymnp1 + emptymnm1 + emptypnm1;

if rn < emptyp1/naybur
    cellsnew(x+1) = 0.5;
    TACage(x+1) = 0;
    
elseif rn < (emptyp1+emptym1)/naybur
    cellsnew(x-1) = 0.5;
    TACage(x-1) = 0;
    
elseif rn < (emptyp1+emptym1+emptypn)/naybur
    cellsnew(x+n) = 0.5;
    TACage(x+n) = 0;
    
elseif rn < (emptyp1+emptym1+emptypn+emptypnp1)/naybur
    cellsnew(x+n+1) = 0.5;
    TACage(x+n+1) = 0;
    
elseif rn < (emptyp1+emptym1+emptypn+emptypnp1+emptymnp1)/naybur
    cellsnew(x+1-n) = 0.5;
    TACage(x+1-n) = 0;
    
elseif rn < (emptyp1+emptym1+emptypn+emptypnp1+emptymnp1+emptymnm1)/naybur
    cellsnew(x-n-1) = 0.5;
    TACage(x-n-1) = 0;
    
elseif rn < (emptyp1+emptym1+emptypn+emptypnp1+emptymnp1+emptymnm1+emptypnm1)/naybur
    cellsnew(x+n-1) = 0.5;
    TACage(x+n-1) = 0;
    
else cellsnew(x-n) = 0.5;
    TACage(x-n) = 0;
    
end

end