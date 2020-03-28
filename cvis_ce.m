function cvis_ce()
    
    format long;
    rng('default');
    
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
    nsamples = 10000;
    
    % limit state function
    g = @(x) Q1(x);
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init,mu,std);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    
%     gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
%     q = @(y) pdf(gm,y);
%     nx = @(x) normpdf(x,mu,std);
%     ny = @(y) normpdf(y,mu,std);
%     probx = integral(nx,l,Inf);
%     proby = integral(ny,l1,Inf);
%     qx = @(x) ((Q(x)<0).*nx(x))./probx;
%     qy = @(y) ((Q1(y)<0).*ny(y))./proby;
%     figure(1)
%     hold on
%     X = linspace(0,5,10000)';
% %     ss = random(gm,5000);
%     plot(X,qx(X),'-',X,qy(X),'--',X,q(X))
% %     histogram(ss,100,'Normalization','pdf');
%     legend('qx','qy','q')
%     hold off
    
%     X = linspace(0,5,1000)';
%     N = size(X,1);
%     h_pre = zeros(N,size(Pi_hat,1));
%     for q = 1:size(Pi_hat,1)
%         h_pre(:,q) = Pi_hat(q)*mvnpdf(X,mu_hat(q,:),Si_hat(:,:,q));
%     end
%     h = sum(h_pre,2);
%     
%     k = (Q(X)<0).*mvnpdf(X,mu,std)/prob;
%     
%     figure(1)
%     hold on
%     plot(X,h,'-o',X,k,'--*')
%     legend('h','k')
%     hold off
%     
%     samples = GM_sample(mu_hat,Si_hat,Pi_hat,100);
%       
%     subplot(1,2,1)
%     plot(X,h,'-o',X,k,'--*')
%     legend('h','k')
%     
%     subplot(1,2,2)
%     hist(samples,100)

    for i = 1:length(a)
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
    for i = 1:length(a)
        e1(i,:) = wQs(i,:)+a(i)*(wQ1s(i,:)-prob1);
        corrcoef(wQs(i,:),wQ1s(i,:))
    end
      
    prob
    prob1
    
    m1 = mean(e1,2);
    m2 = mean(e2,2);
    display('prob./m1')
    prob./m1
    display('prob./m2')
    prob./m2
    
    v1 = var(e1,0,2)
    v2 = var(e2,0,2)
    display('v./v1')
    v./v1
    display('v./v2')
    v./v2
    display('v1./v2')
    v1./v2
    
%     t1 = var(wQs,0,2)
%     t2 = var(wQ1s,0,2)
%     t3 = cov(wQs(1,:),wQ1s(1,:))
    
    figure(1)
    hold on
    plot(a,v1,'-o',a,v2,'--*')
    legend('v_0','v_1')
    xlabel('alpha')
    hold off
    
    figure(2)
    hold on
    plot(a,log(v1),'-o',a,log(v2),'--*')
    legend('ln(v_0)','ln(v_1)')
    xlabel('alpha')
    hold off
    
    figure(3)
    hold on
    plot(a,v1./v2,'-o',a,ones(1,length(a)),'--')
    legend('v_0/v_1')
    xlabel('alpha')
    hold off
    
    astar(1:length(a)) = 0;
    parfor i = 1:length(a)
        cov01 = cov(wQs(i,:),wQ1s(i,:));
        astar(i) = -cov01(1,2)/v2(i);
    end
    astar
    
end