function cvis_ce()
    
    format long;
    % rng('default');
                  
    nelx = 60;
    nely = 20;
    dof = 2*(nely+1)*nelx+2;
    
%     force = 0.01;
%     U = FE(nelx,nely,dof,force);
%     U(dof) 
%     30x10 1.207661533265407
%     60x20 1.228018816903080
    
%     nsamples = 1000000;
%     a = 0.1;
%     b = 0;
%     forces = a+(b-a)*rand(nsamples,1);
%     Ux(1:nsamples) = 0;
%     Uy(1:nsamples) = 0;
%     
%     parfor i = 1:nsamples
%         U = FE(nelx,nely,KE,dof,forces(i));
%         Ux(i) = U(dof-1);
%         Uy(i) = U(dof);
%     end
%     
%     mean(Ux)
%     my = mean(Uy)
%     
%     l = 2*my;
%     mean(l-Uy<0)
    
%     x = linspace(0,nelx,nelx+1);
%     y = linspace(0,nely,nely+1);
%     [X,Y] = meshgrid(x,y);
%     plot (X(:)+Ux,Y(:)+Uy,'o','markersize',4,'markerfacecolor','black');
    
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