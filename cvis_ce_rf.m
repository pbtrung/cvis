function cvis_ce_rf()
    
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
    
    [Psi,lambda,PsiE1] = KL(lx,ly,Lx,Ly,nelx,nely,2);
    outer = 1000;
    inner = 1000;
    EQ_EQ1(lx,ly,nelx,nely,nelx1,nely1,outer,inner,Psi,lambda,PsiE1);

%     neig = length(lambda);
%     L = diag(sqrt(lambda(1:neig)));
%     P = reshape(permute(Psi,[2 1 3]),nelx*nely,neig);
% 
%     l0 = 84.104;
%     l1 = 118.351;
%     EQ0 = 0.001237;
%     VQ0 = 1.194025025025030e-06;
%     EQ1 = 0.020343;
%     
%     dof0 = 2*(nely+1)*nelx+2;
%     dof1 = 2*(nely1+1)*nelx1+2;
%     force = 1;
%     
%     Q0 = @(x) l0-FEQ0(lx,ly,nelx,nely,dof0,force,x,P,L)';
%     Q1 = @(x) l1-FEQ1(lx,ly,nelx1,nely1,dof1,force,x)';
%     
%     a = linspace(-1.5,0.5,33);
%     ansamples = 10000;
%     wQ0s(1:ansamples) = 0;
%     wQ1s(1:ansamples) = 0;
%     
%     % definition of the random variables
%     d      = 1;
%     lower = 1;
%     upper = 2;
%     pi_pdf = repmat(ERADist('uniform','PAR',[lower upper]),d,1);
%     
%     % CE method
%     N      = 5000;    % total number of samples for each level
%     p      = 0.1;     % quantile value to select samples for parameter update
%     k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
%     nsamples = 10000;
%     
%     % limit state function
%     g = @(x) Q1(x);
%     [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init);
%     gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
%     qce = @(x) pdf(gm,x);
%     
%     for j = 1:ansamples
%         samples = random(gm,nsamples);
% 
%         Q0s = Q0(samples(:))<0;
%         Q1s = Q1(samples(:))<0;
% 
%         wQ0s(j) = mean(w.*Q0s);
%         wQ1s(j) = mean(w.*Q1s);
%     end
    
end

function EQ_EQ1(lx,ly,nelx,nely,nelx1,nely1,outer,inner,Psi,lambda,PsiE1)

    force = 1;
%     Uy(1:outer,1:inner) = 0;
%     dof = 2*(nely+1)*nelx+2;
    
    neig = length(lambda);
    L = diag(sqrt(lambda(1:neig)));
    P = reshape(permute(Psi,[2 1 3]),nelx*nely,neig);
%     E(1:nely,1:nelx,1:inner) = 0;
    a = 1;
    b = 2;
        
%     for i = 1:outer
%         for l = 1:inner
%             Z = randn(neig,nelx*nely);
%             LZ = L*Z;
%             PZ = dot(P',LZ);
%             X = reshape(PZ,nelx,nely)';
%             E(:,:,l) = a+(b-a)*normcdf(X);
%         end
%         parfor l = 1:inner
%             fprintf('outer: %d, inner: %d\n',i,l);
%             U = FErf(lx,ly,nelx,nely,dof,force,E(:,:,l));
%             Uy(i,l) = U(dof);
%         end
%     end
%     writematrix(Uy,'Uy.txt');
    
    Uy(1:outer,1:inner) = 0;
    dof = 2*(nely1+1)*nelx1+2;
    for i = 1:outer
        tic
        parfor j = 1:inner
            U = FErf1(lx,ly,nelx1,nely1,dof,force,E1(j));
            Uy(i,j) = U(dof);
        end
        toc
    end
    writematrix(Uy,'Uy1.txt');
    
end