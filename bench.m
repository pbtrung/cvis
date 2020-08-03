function bench()
    
    format long;
    rng('default');
    
    lx = 0.6;
    ly = 0.2;
    nelx = 60;
    nely = 20;
    nelx1 = 30;
    nely1 = 10;
    Lx = 60;
    Ly = 20;
    lower = 1;
    upper = 2;
    
    [Psi,lambda,PsiE1] = KL(lx,ly,Lx,Ly,nelx,nely,2);
    outer = 1;
    inner = 100;
    EQ_EQ1(lx,ly,nelx,nely,nelx1,nely1,outer,inner,Psi,lambda,PsiE1,lower,upper);

end

function EQ_EQ1(lx,ly,nelx,nely,nelx1,nely1,outer,inner,Psi,lambda,PsiE1,a,b)

    t0 = 0;
    t1 = 0;
    force = 1;
    Uy(1:outer,1:inner) = 0;
    dof = 2*(nely+1)*nelx+2;
    
    neig = length(lambda);
    L = diag(sqrt(lambda(1:neig)));
    P = reshape(permute(Psi,[2 1 3]),nelx*nely,neig);
    E(1:nely,1:nelx,1:inner) = 0;
        
    for i = 1:outer
        for j = 1:inner
            Z = randn(neig,1);
            LZ = L*Z;
            PZ = dot(P',repmat(LZ,1,nelx*nely));
            X = reshape(PZ,nelx,nely)';
            E(:,:,j) = a+(b-a)*normcdf(X);
        end
        tic
        for j = 1:inner
%             fprintf('Uy, outer: %d, inner: %d\n',i,j);
            U = FErf(lx,ly,nelx,nely,dof,force,E(:,:,j));
            Uy(i,j) = U(dof);
        end
        t0 = toc
    end
    
    Uy(1:outer,1:inner) = 0;
    dof = 2*(nely1+1)*nelx1+2;
    for i = 1:outer
        E1(1:inner) = 0;
        for j = 1:inner
            Z = randn(neig,1);
            X = PsiE1*(L*Z);
            E1(j) = a+(b-a)*normcdf(X);
        end
        tic
        for j = 1:inner
%             fprintf('Uy1, outer: %d, inner: %d\n',i,j);
            U = FErf1(lx,ly,nelx1,nely1,dof,force,E1(j));
            Uy(i,j) = U(dof);
        end
        t1 = toc
    end
    
    t1/t0
    
end