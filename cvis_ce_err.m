function cvis_ce_err()
    
    format long;
    % rng('default');
    
    mu = 0;
    std = 1;
    b = 3;
    b1 = 2;
    l = mu+b*std;
    l1 = mu+b1*std;
    
    Q = @(x) l-x;
    Q1 = @(y) l1-y;
    
%     vs = 100000;
%     mc(1:vs) = 0;
%     parfor i = 1:vs
%         s = mvnrnd(mu,std,1000000);
%         mc(i) = mean(Q(s(:))<0);
%     end
%     v = var(mc)
%     prob = 1-normcdf(b)
%     prob1 = 1-normcdf(b1)

    v = 1.352806663625048e-09;
    prob = 0.001349898031630;
    prob1 = 0.022750131948179;
    
    % definition of the random variables
    d      = 1;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);
    
    % CE method
    N = [1000 2000 5000 8000];
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
    nsamples = 10000;
    % limit state function
    g = @(x) Q1(x);
    
    ansamples = 100000;
    a(1:length(N)) = 0;
    kldiv(1:length(N)) = 0;
    e1(1:length(N),1:ansamples) = 0;
    e2(1:length(N),1:ansamples) = 0;
    wQs(1:length(N),1:ansamples) = 0;
    wQ1s(1:length(N),1:ansamples) = 0;
    
    n1 = @(y) normpdf(y,mu,std);
    m = @(y) (Q1(y)<0).*n1(y);
    qopt = @(y) m(y)./prob1;
    umin = 0;
    umax = 7;
    nrs = 10000000;
    
    for i = 1:length(N)
        [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N(i),p,g,pi_pdf,k_init,mu,std);
        gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
        qce = @(y) pdf(gm,y);
        
%         y = GM_sample(mu_hat,Si_hat,Pi_hat,10000000);
%         accepted = (Q1(y)<0) ~= 0;
%         s = y(accepted);
        
        vs = 1000;
        mc(1:vs) = 0;
        parfor k = 1:vs
            u = umin+(umax-umin)*rand(nrs,1);
            sample_value = qopt(u);
            max_value = max(sample_value);
            accepted = rand(nrs,1)<(sample_value/max_value);
            s = u(accepted,:);
            mc(k) = mean(log(qopt(s(:)))-log(qce(s(:))));
        end
        kldiv(i) = mean(mc)
        
%         integral(@(y) qopt(y).*(log(qopt(y))-log(qce(y))),-Inf,Inf,'ArrayValued',true,'RelTol',0,'AbsTol',1e-9)
        
        parfor j = 1:ansamples
            samples = random(gm,nsamples);
            q = pdf(gm,samples);
            
            Qs = Q(samples(:))<0;
            Q1s = Q1(samples(:))<0;
            w = mvnpdf(samples,mu,std)./q;

            wQs(i,j) = mean(w.*Qs);
            wQ1s(i,j) = mean(w.*Q1s);
        end
    end
      
    e2 = wQs;
    v2 = var(e2,0,2);
    for i = 1:length(N)
        cov12 = cov(wQs(i,:),wQ1s(i,:));
        a(i) = -cov12(1,2)/v2(i);
        e1(i,:) = wQs(i,:)+a(i)*(wQ1s(i,:)-prob1);
        corrcoef(wQs(i,:),wQ1s(i,:))
    end
    v1 = var(e1,0,2);
    
    prob
    prob1
    v1
    v2
    a
    kldiv
    
    display('prob./m1')
    prob./mean(e1,2)
    display('prob./m2')
    prob./mean(e2,2)
    
    display('v./v1')
    v./v1
    display('v./v2')
    v./v2
    display('v1./v2')
    v1./v2
    
    figure(1)
    hold on
    plot(N,v1./v2,'-o')
    legend('v_0/v_1')
    xlabel('N')
    hold off
    
    figure(2)
    hold on
    plot(N,log(v1./v2),'-o')
    legend('ln(v_0/v_1)')
    xlabel('N')
    hold off
    
    figure(3)
    hold on
    plot(N,kldiv,'-o')
    legend('KL Divergence')
    xlabel('N')
    hold off
    
end