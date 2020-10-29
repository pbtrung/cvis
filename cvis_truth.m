function cvis_truth(path, model)
    
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
    
    [Psi,lambda,~] = KL(lx,ly,Lx,Ly,nelx0,nely0,fl);
    
    ns = 10^6;
    Uy(1:ns) = 0;
    force = 1;
    dof0 = 2*(nely0+1)*nelx0+2;
    dof1 = 2*(nely1+1)*nelx1+2;
    neig = length(lambda);
    L = diag(sqrt(lambda(1:neig)));
    P = reshape(permute(Psi,[2 1 3]),nelx0*nely0,neig);
    Z = randn(neig,ns);
    
    if model == 0
        parfor i = 1:ns
            LZ = L*Z(:,i);
            PZ = dot(P',repmat(LZ,1,nelx0*nely0));
            X = reshape(PZ,nelx0,nely0)';
            E = lower + (upper-lower)*normcdf(X);
            tic
            U = FE(lx,ly,nelx0,nely0,dof0,force,E,model);
            t = toc;
            Uy(i) = U(dof0);
            fprintf('iter: %d, time: %f\n',i,t);
        end
        writematrix(Uy,append(path,'Uy0.txt'));
    elseif model == 1
        parfor i = 1:ns
            X = PsiE1*(L*Z(:,i));
            E1 = lower + (upper-lower)*normcdf(X);
            tic
            U = FE(lx,ly,nelx1,nely1,dof1,force,E1,model);
            t = toc;
            Uy(i) = U(dof1);
            fprintf('iter: %d, time: %f\n',i,t);
        end
        writematrix(Uy,append(path,'Uy1.txt'));
    end
    
end