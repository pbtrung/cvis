function cvis_ce()
    E = 10000;
    poisson = 0.30;
    kapa = 5 / 6;
    L = 1;
    
    nelx = 40;
    nely = 40;
    nelx1 = 10;
    nely1 = 10;
    if nelx ~= nely || mod(nelx,2) ~= 0 || nelx1 ~= nely1 || mod(nelx1,2) ~= 0
        error("Check inputs")
    end
    
    numElems = nelx*nely;
    numElems1 = nelx1*nely1;
    
end

