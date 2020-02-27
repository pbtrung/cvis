function cvis_ce()
    
    format long;
    rng('default');
    
    mu = 0;
    std = 1;
    l = mu+3*std;
    l1 = mu+2*std;
    
    Q = @(x) l-x;
    Q1 = @(y) l1-y;
    
    s = mvnrnd(mu,std,10000000);
    prob = mean(Q(s(:))<0);
    v = var(Q(s(:))<0);
    prob1 = mean(Q1(s(:))<0);
    
    a = linspace(-4,4,33);
    ansamples = 100000;
    e1(1:length(a),ansamples) = 0;
    e2(1:length(a),ansamples) = 0;
    wQs(1:length(a),ansamples) = 0;
    wQ1s(1:length(a),ansamples) = 0;
    
    % definition of the random variables
    d      = 1;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);

    % CE method
    N      = 5000;    % total number of samples for each level
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
    nsamples = 100000;
    
    % limit state function
    g = @(x) Q1(x);
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init,mu,std);
    
    for i = 1:length(a)
        parfor j = 1:ansamples
            samples = GM_sample(mu_hat,Si_hat,Pi_hat,nsamples);
            q = q_calc(samples,mu_hat,Si_hat,Pi_hat);

            Qs = Q(samples(:))<0;
            Q1s = Q1(samples(:))<0;
            w = mvnpdf(samples,mu,std)./q;

            wQs(i,j) = mean(w.*Qs);
            wQ1s(i,j) = mean(w.*Q1s);
        end
    end
    
    e2 = wQs;
    for i = 1:length(a)
        e1(i,:) = wQs(i,:)+a(i)*(wQ1s(i,:)-prob1);
        corrcoef(wQs(i,:),wQ1s(i,:))
    end
    
    prob
    prob1
    
    m1 = mean(e1,2);
    m2 = mean(e2,2);
    prob./m1
    prob./m2
    
    v1 = var(e1,0,2);
    v2 = var(e2,0,2);
    v./v1
    v./v2
    v1./v2
    
    figure(1)
    hold on
    plot(a,v1,'-o',a,v2,'--*')
    legend('v1','v2')
    xlabel('alpha')
    hold off
    
    figure(2)
    hold on
    plot(a,log(v1),'-o',a,log(v2),'--*')
    legend('log(v1)','log(v2)')
    xlabel('alpha')
    hold off
    
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
