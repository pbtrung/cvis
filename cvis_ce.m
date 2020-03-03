function cvis_ce()
    
    format long;
    rng('default');
    
    lx = 60;
    ly = 20;
    nelx = 60;
    nely = 20;
    nelx1 = 30;
    nely1 = 10;
    mu = 0;
    sigma = 1;
%     nsamples = 1000000;
%     [EQ,EQ1,VQ] = EQ_EQ1(lx,ly,nelx,nely,nelx1,nely1,nsamples,mu,sigma)

    EQ = 0.001251;
    EQ1 = 0.031769;
    VQ = 0.001249436248438;
    l = 3.705226467613969e+02;
    l1 = 2.242340722569870e+02;
    
    KE = ElmStiffnessMatrix(lx,ly,nelx,nely);
    KE1 = ElmStiffnessMatrix(lx,ly,nelx1,nely1);
    dof = 2*(nely+1)*nelx+2;
    dof1 = 2*(nely1+1)*nelx1+2;
    
    Qx = @(x) l-FE3(nelx,nely,dof,KE,x);
    Q1y = @(y) l1-FE3(nelx1,nely1,dof1,KE1,y);
    select = @(v,idx) v(idx,:);
    Q = @(x) select(Qx(x),dof)';
    Q1 = @(y) select(Q1y(y),dof1)';
    
    alpha = linspace(-4,4,33);
    ansamples = 100;
    e1(1:length(alpha),ansamples) = 0;
    e2(1:length(alpha),ansamples) = 0;
    wQs(1:length(alpha),ansamples) = 0;
    wQ1s(1:length(alpha),ansamples) = 0;
    
    % definition of the random variables
    d      = 1;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);

    % CE method
    N      = 5000;    % total number of samples for each level
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
    nsamples = 10000;
    % nsamples=10000 needs 110s per core
    
    % limit state function
    g = @(x) Q1(x);
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init,mu,sigma);
    
    tic
    for i = 1:length(alpha)
        parfor j = 1:ansamples
            fprintf('alpha: %f, iter: %d\n',alpha(i),j);
            samples = GM_sample(mu_hat,Si_hat,Pi_hat,nsamples);
            q = q_calc(samples,mu_hat,Si_hat,Pi_hat);

            Qs = Q(samples(:))<0;
            Q1s = Q1(samples(:))<0;
            w = mvnpdf(samples,mu,sigma)./q;

            wQs(i,j) = mean(w.*Qs);
            wQ1s(i,j) = mean(w.*Q1s);
        end
        fprintf('alpha: %f, time: %f\n',alpha(i),toc);
    end
    fprintf('finish: %f\n',toc);
    
    e2 = wQs;
    for i = 1:length(alpha)
        e1(i,:) = wQs(i,:)+alpha(i)*(wQ1s(i,:)-EQ1);
        corrcoef(wQs(i,:),wQ1s(i,:))
    end
    
    EQ
    EQ1
    
    m1 = mean(e1,2);
    m2 = mean(e2,2);
    EQ./m1
    EQ./m2
    
    v1 = var(e1,0,2);
    v2 = var(e2,0,2);
    VQ./v1
    VQ./v2
    v1./v2
    
    figure(1)
    hold on
    plot(alpha,v1,'-o',alpha,v2,'--*')
    legend('v1','v2')
    xlabel('alpha')
    hold off
    
    figure(2)
    hold on
    plot(alpha,log(v1),'-o',alpha,log(v2),'--*')
    legend('log(v1)','log(v2)')
    xlabel('alpha')
    hold off
    
end

function [EQ,EQ1,VQ] = EQ_EQ1(lx,ly,nelx,nely,nelx1,nely1,nsamples,mu,sigma)
    forces = mvnrnd(mu,sigma,nsamples);
    Uy(1:nsamples) = 0;
    
    dof = 2*(nely+1)*nelx+2;
    KE = ElmStiffnessMatrix(lx,ly,nelx,nely);
    parfor i = 1:nsamples
        U = FE2(nelx,nely,dof,KE,forces(i));
        Uy(i) = U(dof);
    end
    m = max(Uy);
    my = mean(Uy);
    l = 0.65*m
    EQ = mean(l-Uy<0);
    VQ = var(l-Uy<0);
    
    Uy(1:nsamples) = 0;
    dof = 2*(nely1+1)*nelx1+2;
    KE = ElmStiffnessMatrix(lx,ly,nelx1,nely1);
    parfor i = 1:nsamples
        U = FE2(nelx1,nely1,dof,KE,forces(i));
        Uy(i) = U(dof);
    end
    m = max(Uy);
    my = mean(Uy);
    l = 0.4*m
    EQ1 = mean(l-Uy<0);
end

function KE = ElmStiffnessMatrix(lx,ly,nelx,nely)
    E = 1.;
    nu = 0.3;
    h = 1;
    a = lx/nelx;
    b = ly/nely;
    c = b/a;
    
    ka = (2/c)*(1-nu);
    kb = (3/2)*(1+nu);
    kc = (3/2)*(1-3*nu);
    kd = (2*c)*(1-nu);
    KE = (E*h)/(12*(1-nu^2))*[ ...
        4*c+ka kb -4*c+ka/2 -kc -2*c-ka/2 -kb 2*c-ka kc
        kb 4/c+kd kc 2/c-kd -kb -2/c-kd/2 -kc -4/c+kd/2
        -4*c+ka/2 kc 4*c+ka -kb 2*c-ka -kc -2*c-ka/2 kb
        -kc 2/c-kd -kb 4/c+kd kc -4/c+kd/2 kb -2/c-kd/2
        -2*c-ka/2 -kb 2*c-ka kc 4*c+ka kb -4*c+ka/2 -kc
        -kb -2/c-kd/2 -kc -4/c+kd/2 kb 4/c+kd kc 2/c-kd
        2*c-ka -kc -2*c-ka/2 kb -4*c+ka/2 kc 4*c+ka -kb
        kc -4/c+kd/2 kb -2/c-kd/2 -kc 2/c-kd -kb 4/c+kd
        ];
end

function X = GM_sample(mu,Si,Pi,N)
    % Algorithm to draw samples from a Gaussian-Mixture (GM) distribution
    %{
    ---------------------------------------------------------------------------
    Input:
    * mu : [npi x d]-array of means of Gaussians in the Mixture
    * Si : [d x d x npi]-array of cov-matrices of Gaussians in the Mixture
    * Pi : [npi]-array of weights of Gaussians in the Mixture (sum(Pi) = 1)
    * N  : number of samples to draw from the GM distribution
    ---------------------------------------------------------------------------
    Output:
    * X  : samples from the GM distribution
    ---------------------------------------------------------------------------
    %}

    if size(mu,1) == 1
        X = mvnrnd(mu,Si,N);
    else
        % Determine number of samples from each distribution

        ind = randsample(size(mu,1),N,true,Pi);
        z = histcounts(ind,[(1:size(mu,1)) size(mu,1)+1]);
        % Generate samples
        X   = zeros(N,size(mu,2));
        ind = 1;
        for i = 1:size(Pi,1)
            np                = z(i);
            X(ind:ind+np-1,:) = mvnrnd(mu(i,:),Si(:,:,i),np);
            ind               = ind+np;
        end
    end
end

function h = q_calc(X,mu,Si,Pi)
    % Basic algorithm to calculate h for the likelihood ratio
    %{
    ---------------------------------------------------------------------------
    Input:
    * X  : input samples
    * mu : [npi x d]-array of means of Gaussians in the Mixture
    * Si : [d x d x npi]-array of cov-matrices of Gaussians in the Mixture
    * Pi : [npi]-array of weights of Gaussians in the Mixture (sum(Pi) = 1)
    ---------------------------------------------------------------------------
    Output:
    * h  : parameters h (IS density)
    ---------------------------------------------------------------------------
    %}

    N = size(X,1);
    if size(Pi,1) == 1
        h = mvnpdf(X,mu,Si);
    else
        h_pre = zeros(N,size(Pi,1));
        for q = 1:size(Pi,1)
            h_pre(:,q) = Pi(q)*mvnpdf(X,mu(q,:),Si(:,:,q));
        end
        h = sum(h_pre,2);
    end
end