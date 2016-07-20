function [cellsnew,TACage] = PlaceTACDaughter8Speedup(cellsnew,x,rn,n,TACage)
% this function places a daughter of the other type into a random empty neighbouring location

[emptyp1,emptym1,emptypn,emptymn,emptypnp1,emptymnp1,emptymnm1,emptypnm1] = SpaceCheck8(x,cellsnew,n); %check for space
naybur = emptyp1 + emptym1 + emptypn + emptymn + emptypnp1 + emptymnp1 + emptymnm1 + emptypnm1;

if rn < emptyp1/naybur %was if, now elseif
    cellsnew(x+1) = 1;
    TACage(x+1) = TACage(x)+1;
elseif rn < (emptyp1+emptym1)/naybur
    cellsnew(x-1) = 1;
    TACage(x-1) = TACage(x)+1;
elseif rn < (emptyp1+emptym1+emptypn)/naybur
    cellsnew(x+n) = 1;
    TACage(x+n) = TACage(x)+1;
elseif rn < (emptyp1+emptym1+emptypn+emptypnp1)/naybur
    cellsnew(x+n+1) = 1;
    TACage(x+n+1) = TACage(x)+1;
elseif rn < (emptyp1+emptym1+emptypn+emptypnp1+emptymnp1)/naybur
    cellsnew(x+1-n) = 1;
    TACage(x+1-n) = TACage(x)+1;
elseif rn < (emptyp1+emptym1+emptypn+emptypnp1+emptymnp1+emptymnm1)/naybur
    cellsnew(x-n-1) = 1;
    TACage(x-n-1) = TACage(x)+1;
elseif rn < (emptyp1+emptym1+emptypn+emptypnp1+emptymnp1+emptymnm1+emptypnm1)/naybur
    cellsnew(x+n-1) = 1;
    TACage(x+n-1) = TACage(x)+1;
else cellsnew(x-n) = 1;
    TACage(x-n) = TACage(x)+1;
end

