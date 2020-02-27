function cvis_ce()
    
    format long;
    % rng('default');
    
    lx = 60;
    ly = 20;
    nelx = 60;
    nely = 20;
    dof = 2*(nely+1)*nelx+2;
    KE = ElmStiffnessMatrix(lx,ly,nelx,nely);
       
%     force = 0.05;
%     U = FE2(nelx,nely,dof,KE,force);
%     U(dof)
    
%     U = FE1(nelx,nely,dof,force);
%     U(dof)
        
    nsamples = 100000;
    a = 0.1;
    b = 0;
    forces = a+(b-a)*rand(nsamples,1);
    Ux(1:nsamples) = 0;
    Uy(1:nsamples) = 0;
    
    tic
    parfor i = 1:nsamples
        U = FE2(nelx,nely,dof,KE,forces(i));
        Ux(i) = U(dof-1);
        Uy(i) = U(dof);
    end
    toc
    
    mean(Ux)
    my = mean(Uy)
    
    l = 1.995*my;
    mean(l-Uy<0)
    
%     x = linspace(0,nelx,nelx+1);
%     y = linspace(0,nely,nely+1);
%     [X,Y] = meshgrid(x,y);
%     plot (X(:)+Ux,Y(:)+Uy,'o','markersize',4,'markerfacecolor','black');
    
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