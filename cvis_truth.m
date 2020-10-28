function cvis_truth(path)
    
    format long;
    rng('default');
    
    lx = 0.6;
    ly = 0.2;
    nelx0 = 60;
    nely0 = 20;
    Lx = 60;
    Ly = 20;
    lower = 1;
    upper = 2;
    fl = 2;
    
    [Psi,lambda,~] = KL(lx,ly,Lx,Ly,nelx0,nely0,fl);
    
    ns = 10^6;
    Uy(1:ns) = 0;
    force = 1;
    dof0 = 2*(nely0+1)*nelx0+2;
    neig = length(lambda);
    L = diag(sqrt(lambda(1:neig)));
    P = reshape(permute(Psi,[2 1 3]),nelx0*nely0,neig);
    Z = randn(neig,ns);
    
    parfor i = 1:ns
        LZ = L*Z(:,i);
        PZ = dot(P',repmat(LZ,1,nelx0*nely0));
        X = reshape(PZ,nelx0,nely0)';
        E = lower + (upper-lower)*normcdf(X);
        U = FE(lx,ly,nelx0,nely0,dof0,force,E,0);
        Uy(i) = U(dof0);
    end
    
    writematrix(Uy,append(path,'Uy.txt'));
    
end