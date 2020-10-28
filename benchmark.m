function benchmark()
    
    format long;
    rng('default');
    
    lx = 0.6;
    ly = 0.2;
    nelx0 = 60;
    nely0 = 20;
    nelx1 = 30;
    nely1 = 10;
    Lx = 60;
    Ly = 20;
    lower = 1;
    upper = 2;
    fl = 2;
    
    [Psi,lambda,PsiE1] = KL(lx,ly,Lx,Ly,nelx0,nely0,fl);
    
    ns = 1000;
    force = 1;
    dof0 = 2*(nely0+1)*nelx0+2;
    dof1 = 2*(nely1+1)*nelx1+2;
    neig = length(lambda);
    L = diag(sqrt(lambda(1:neig)));
    P = reshape(permute(Psi,[2 1 3]),nelx0*nely0,neig);
    E(1:nely0,1:nelx0,1:ns) = 0;
    
    t0(1:ns) = 0;
    t1(1:ns) = 0;
    
    for i = 1:ns
        Z = randn(neig,1);
        LZ = L*Z;
        PZ = dot(P',repmat(LZ,1,nelx0*nely0));
        X = reshape(PZ,nelx0,nely0)';
        E(:,:,i) = lower + (upper-lower)*normcdf(X);
        tic
        FE(lx,ly,nelx0,nely0,dof0,force,E(:,:,i),0);
        t0(i) = toc;
    end
    
    for i = 1:ns
        Z = randn(neig,1);
        X = PsiE1*(L*Z);
        E1 = lower + (upper-lower)*normcdf(X);
        tic
        FE(lx,ly,nelx1,nely1,dof1,force,E1,1);
        t1(i) = toc;
    end
    
    disp(mean(t0)/mean(t1));
    
end