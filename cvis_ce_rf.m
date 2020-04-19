function cvis_ce_rf()
    
    format long;
    rng('default');
    
    lx = 0.6;
    ly = 0.2;
    nelx = 60;
    nely = 20;
    nelx1 = 30;
    nely1 = 10;
    Lx = 0.6;
    Ly = 0.2;
    
    [Psi,lambda] = KL(lx,ly,Lx,Ly,nelx,nely);
    outer = 1000;
    inner = 1000;
    [EQ,EQ1,VQ] = EQ_EQ1(lx,ly,nelx,nely,nelx1,nely1,outer,inner,Psi,lambda);
    
end

function [EQ,EQ1,VQ] = EQ_EQ1(lx,ly,nelx,nely,nelx1,nely1,outer,inner,Psi,lambda)

    force = 1;
    Uy(1:outer,1:inner) = 0;
    dof = 2*(nely+1)*nelx+2;
    
    neig = length(lambda);
    L = diag(sqrt(lambda(1:neig)));
    E(1:nely,1:1:nelx,1:inner) = 0;    
    a = 1;
    b = 2;
    
    for i = 1:outer
        for l = 1:inner
            Z = randn(neig,nelx*nely);
            LZ = L*Z;
            for j = 1:nelx
                for k = 1:nely
                    X = squeeze(Psi(k,j,:))'*LZ(:,(k-1)*nelx+j);
                    E(k,j,l) = a+(b-a)*normcdf(X);
                end
            end
        end
        parfor l = 1:inner
            U = FErf(lx,ly,nelx,nely,dof,force,E(:,:,l));
            Uy(i,l) = U(dof);
        end
    end
    
    m = max(Uy(:));
    my = mean(mean(Uy,2));
    l = 0.65*m
    EQ = mean(mean(l-Uy<0,2));
    VQ = var(mean(l-Uy<0,2));
    
    Uy(1:outer,1:inner) = 0;
    dof = 2*(nely1+1)*nelx1+2;
    E1(1:nely1,1:1:nelx1,1:inner) = 0;
    for i = 1:outer
        for l = 1:inner
            Z = randn(neig,nelx1*nely1);
            LZ = L*Z;
            for j = 1:nelx1
                for k = 1:nely1
                    X = squeeze(Psi(k,j,:))'*LZ(:,(k-1)*nelx1+j);
                    E1(k,j,l) = a+(b-a)*normcdf(X);
                end
            end
        end
        parfor l = 1:inner
            U = FErf(lx,ly,nelx1,nely1,dof,force,E1(:,:,l));
            Uy(i,l) = U(dof);
        end
    end
    m = max(Uy(:));
    my = mean(mean(Uy,2));
    l = 0.4*m
    EQ1 = mean(mean(l-Uy<0,2));
    
end