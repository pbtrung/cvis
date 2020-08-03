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
%     outer = 1000;
%     inner = 1000;
%     EQ_EQ1(lx,ly,nelx,nely,nelx1,nely1,outer,inner,Psi,lambda,PsiE1,lower,upper);

    neig = length(lambda);
    L = diag(sqrt(lambda(1:neig)));
    P = reshape(permute(Psi,[2 1 3]),nelx*nely,neig);

    l0 = 118.923;
    l1 = 108.510;
    
    dof0 = 2*(nely+1)*nelx+2;
    dof1 = 2*(nely1+1)*nelx1+2;
    force = 1;
    
    Q0 = @(x) l0-FEQ0(lx,ly,nelx,nely,dof0,force,x,P,L,lower,upper)';
    Q1 = @(x) l1-FEQ1(lx,ly,nelx1,nely1,dof1,force,x,PsiE1,L,lower,upper)';
        
    K = 100;
    nsamples = 100;
    
    % definition of the random variables
    d      = neig;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);
    
    % CE method
%     N      = 5000;    % total number of samples for each level
%     p      = 0.1;     % quantile value to select samples for parameter update
%     k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
%      
%     % limit state function
%     g = @(x) Q1(x);
%     [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init);
%     gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    
    for j = 1:2
        size(random(gm,nsamples))
%         Q0s = Q0(samples')<0;
%         Q1s = Q1(samples')<0;
    end

end