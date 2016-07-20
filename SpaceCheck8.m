function [emptyp1,emptym1,emptypn,emptymn,emptypnp1,emptymnp1,emptymnm1,emptypnm1] = SpaceCheck8(x,cellarray,n)

%check neighborhood for empty space

%corners
if x == 1 %top left corner
    
            emptym1 = 0;
            emptymn = 0;
            emptymnp1 = 0;
            emptymnm1 = 0;
            emptypnm1 = 0;
    
    if cellarray(x+1) == 0
        emptyp1 = 1;
    else emptyp1 = 0;
    end
    
    if cellarray(x+n) == 0
        emptypn = 1;
    else emptypn = 0;
    end
    
    if cellarray(x+n+1) == 0
        emptypnp1 = 1;
    else emptypnp1 = 0;
    end
    
elseif x == n % top right corner
    
    emptymn = 0;
    emptyp1 = 0;
    emptypnp1 = 0;
    emptymnp1 = 0;
    emptymnm1 = 0;

    if cellarray(x-1) == 0
        emptym1 = 1;
    else emptym1 = 0;
    end
    
    if cellarray(x+n) == 0
        emptypn = 1;
    else emptypn = 0;
    end
    
    if cellarray(x+n-1) == 0
        emptypnm1 = 1;
    else emptypnm1 = 0;
    end
    
elseif x == n*n % bottom right corner
    
    emptyp1 = 0;
    emptypn = 0;
    emptypnp1 = 0;
    emptymnp1 = 0;
    emptypnm1 = 0;
    
    if cellarray(x-1) == 0
        emptym1 = 1;
    else emptym1 = 0;
    end
    
    if cellarray(x-n) == 0
        emptymn = 1;
    else emptymn =0;
    end
    
    if cellarray(x-n-1) == 0
        emptymnm1 = 1;
    else emptymnm1 = 0;
    end
    
elseif x == ((n-1)*n + 1) % bottom left corner
    
    emptym1 = 0;
    emptypn = 0;
    emptypnp1 = 0;
    emptymnm1 = 0;
    emptypnm1 = 0;
    
    if cellarray(x+1) == 0
        emptyp1 = 1;
    else emptyp1 = 0;
    end
    
    if cellarray(x-n) == 0
        emptymn = 1;
    else emptymn =0;
    end
    
    if cellarray(x-n+1) == 0
        emptymnp1 = 1;
    else emptymnp1 = 0;
    end
    
elseif x <= n %if you are in the top row, dont' look up!
    
    emptymn = 0;
    emptymnp1 = 0;
    emptymnm1 = 0;
    
    if cellarray(x+1) == 0
        emptyp1 = 1;
    else emptyp1 = 0;
    end
    
    if cellarray(x-1) == 0
        emptym1 = 1;
    else emptym1 = 0;
    end
    
    if cellarray(x+n) == 0
        emptypn = 1;
    else emptypn = 0;
    end
    
    if cellarray(x+n+1) == 0
        emptypnp1 = 1;
    else emptypnp1 = 0;
    end
    
    if cellarray(x+n-1) == 0
        emptypnm1 = 1;
    else emptypnm1 = 0;
    end
    
elseif x > n*(n-1) % if you are in the bottom row, don't look down
    
    emptypn = 0;
    emptypnp1 = 0;
    emptypnm1 = 0;
    
    if cellarray(x+1) == 0
        emptyp1 = 1;
    else emptyp1 = 0;
    end
    
    if cellarray(x-1) == 0
        emptym1 = 1;
    else emptym1 = 0;
    end

    if cellarray(x-n) == 0
        emptymn = 1;
    else emptymn =0;
    end

    if cellarray(x-n+1) == 0
        emptymnp1 = 1;
    else emptymnp1 = 0;
    end
    
    if cellarray(x-n-1) == 0
        emptymnm1 = 1;
    else emptymnm1 = 0;
    end

elseif mod(x-1,n) == 0 % if you are in the left hand row, don't look left
    
    emptym1 = 0;
    emptymnm1 = 0;
    emptypnm1 = 0;
    
    if cellarray(x+1) == 0
        emptyp1 = 1;
    else emptyp1 = 0;
    end
    
    if cellarray(x+n) == 0
        emptypn = 1;
    else emptypn =0;
    end
    
    if cellarray(x-n) == 0
        emptymn = 1;
    else emptymn =0;
    end
    
    if cellarray(x+n+1) == 0
        emptypnp1 = 1;
    else emptypnp1 = 0;
    end
    
    if cellarray(x-n+1) == 0
        emptymnp1 = 1;
    else emptymnp1 = 0;
    end

elseif mod(x,n) == 0 % if you are in the right hand row, don't look right
    
    emptyp1 = 0;
    emptypnp1 = 0;
    emptymnp1 = 0;
 
    if cellarray(x-1) == 0
        emptym1 = 1;
    else emptym1 = 0;
    end
    
    if cellarray(x+n) == 0
        emptypn = 1;
    else emptypn = 0;
    end
    
    if cellarray(x-n) == 0
        emptymn = 1;
    else emptymn = 0;
    end
  
    if cellarray(x-n-1) == 0
        emptymnm1 = 1;
    else emptymnm1 = 0;
    end
    
    if cellarray(x+n-1) == 0
        emptypnm1 = 1;
    else emptypnm1 = 0;
    end
    
else % not on an edge or in a corner
    
    if cellarray(x+1) == 0
        emptyp1 = 1;
    else emptyp1 = 0;
    end
    
    if cellarray(x-1) == 0
        emptym1 = 1;
    else emptym1 = 0;
    end
    
    if cellarray(x+n) == 0
        emptypn = 1;
    else emptypn = 0;
    end
    
    if cellarray(x-n) == 0
        emptymn = 1;
    else emptymn = 0;
    end
    
    if cellarray(x+n+1) == 0
        emptypnp1 = 1;
    else emptypnp1 = 0;
    end
    
    if cellarray(x-n+1) == 0
        emptymnp1 = 1;
    else emptymnp1 =0;
    end
    
    if cellarray(x-n-1) == 0
        emptymnm1 = 1;
    else emptymnm1 = 0;
    end
    
    if cellarray(x+n-1) == 0
        emptypnm1 = 1;
    else emptypnm1 = 0;
    end
    
end