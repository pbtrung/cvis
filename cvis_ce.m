function cvis_ce()    

    format long;
    rng('default');

    nelx = 30;
    nely = 30;
    nelx1 = 10;
    nely1 = 10;
    if nelx ~= nely || mod(nelx,2) ~= 0 || nelx1 ~= nely1 || mod(nelx1,2) ~= 0
        error("Check inputs")
    end
        
%     l0 = 
    l1 = 0.010287;
%     EQ0 = 
%     VQ0 = 
    EQ1 = 0.019249;
    
    % x(nsamples,8)
    Q0 = @(x) l0-FE_plate(nelx,nely,x);
    Q1 = @(x) l1-FE_plate(nelx1,nely1,x);
    
    a = linspace(-1.5,0.5,33);
    ansamples = 100;
    wQ0s(1:ansamples) = 0;
    wQ1s(1:ansamples) = 0;
    
    % definition of the random variables
    d = 8;
    lower = 1;
    upper = 2;
    pi_pdf = repmat(ERADist('uniform','PAR',[lower upper]),d,1);
    if any(strcmp('Marginals',fieldnames(pi_pdf))) == 1
        % use Nataf transform (dependence)
        u2x = @(u) pi_pdf.U2X(u);
    else
        % use distribution information for the transformation (independence)
        u2x = @(u) pi_pdf(1).icdf(normcdf(u));
    end
    
    % CE method
    N      = 5000;    % total number of samples for each level
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
    nsamples = 100;
    
    % limit state function
    g = @(x) Q1(x');
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    qce = @(x) pdf(gm,x);
    
    for j = 1:ansamples
        fprintf('ansamples: %d\n',j);
        tic
        samples = random(gm,nsamples);
        
        w = mvnpdf(samples,zeros(1,d),eye(d))./qce(samples);
        samples = u2x(samples);
        Q0s = Q0(samples)<0;
        Q1s = Q1(samples)<0;

        wQ0s(j) = mean(w.*Q0s);
        wQ1s(j) = mean(w.*Q1s);
        toc
    end
    writematrix(wQ0s,'wQ0s.txt');
    writematrix(wQ1s,'wQ1s.txt');
    
%     m0 = mean(wQ0s);
%     m1 = mean(wQ1s);
%     me = m0+a*(m1-EQ1);
      
%     EQ0
    EQ1
    
%     display('EQ0./m0')
%     EQ0./m0
%     display('EQ0./me')
%     EQ0./me'
    
    covar = cov(wQ0s,wQ1s);
    v0 = covar(1,1);
    v1 = covar(2,2);
    ve = v0+a.^2*v1+2*a*covar(1,2);
      
%     display('VQ0/v0')
%     VQ0/v0
%     display('VQ0./ve')
%     VQ0./ve'
    display('ve./v0')
    ve'./v0
    
    figure(1)
    hold on
    plot(a,ve,'-o',a,v0*ones(1,length(a)),'--*')
    legend('v_e','v_0')
    xlabel('alpha')
    hold off
    
    figure(2)
    hold on
    plot(a,ve./v0,'-o',a,ones(1,length(a)),'--*')
    legend('v_e/v_0')
    xlabel('alpha')
    hold off
end