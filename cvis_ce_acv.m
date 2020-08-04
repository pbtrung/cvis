function cvis_ce_acv(path)
    
    format long;
    rng('default');
    
    mu = 0;
    std = 1;
    w0 = 3;
    w1 = 2.8;
    l0 = mu+w0*std;
    l1 = mu+w1*std;
    
    Q0 = @(x) l0-x;
    Q1 = @(x) l1-x;

    v = 1.352806663625048e-09;
    prob0 = 0.001349898031630;
    prob1 = 1-normcdf(l1);
    
    a = linspace(-1.5,0.5,33);
    K = 10000;
    wQ0s(1:K) = 0;
    wQ0s_cv(1:K) = 0;
    wQ1s_cv(1:K) = 0;
    wQ0s_acv1(1:K) = 0;
    wQ1s_acv1(1:K) = 0;
    mu_acv1(1:K) = 0;
    wQ0s_acv2(1:K) = 0;
    wQ1s_acv2(1:K) = 0;
    mu_acv2(1:K) = 0;
    
    % definition of the random variables
    d      = 1;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);
    
    % CE method
    N      = 3000;    % total number of samples for each level
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
    
    ns_mfis = 1000;
    ns_Q0_cv = 910;
    ns_Q1_cv = 910;
    
    ns_Q0_acv1 = 819;
    ns_Q1_acv1 = 819;
    ns_mu_acv1 = 1000;
    
    ns_Q0_acv2 = 819;
    ns_Q1_acv2 = 819;
    ns_mu_acv2 = 1000;
        
    % limit state function
    g = @(x) Q1(x);
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    qce = @(x) pdf(gm,x);

    for j = 1:K
        ss = random(gm,ns_mfis);
        Q0s = Q0(ss(:))<0;
        w = mvnpdf(ss,mu,std)./qce(ss);
        wQ0s(j) = mean(w.*Q0s);
        
        wQ0s_cv(j) = mean(w(1:ns_Q0_cv).*Q0s(1:ns_Q0_cv));
        Q1s = Q1(ss(1:ns_Q1_cv))<0;
        wQ1s_cv(j) = mean(w(1:ns_Q1_cv).*Q1s(1:ns_Q1_cv));
        
        wQ0s_acv1(j) = mean(w(1:ns_Q0_acv1).*Q0s(1:ns_Q0_acv1));
        Q1s_acv1 = Q1(ss(1:ns_Q1_acv1))<0;
        wQ1s_acv1(j) = mean(w(1:ns_Q1_acv1).*Q1s_acv1(1:ns_Q1_acv1));
        ss_mu_acv1 = random(gm,ns_mu_acv1);
        Q1s_mu_acv1 = Q1(ss_mu_acv1)<0;
        w_mu_acv1 = mvnpdf(ss_mu_acv1,mu,std)./qce(ss_mu_acv1);
        mu_acv1(j) = mean(w_mu_acv1(1:ns_mu_acv1).*Q1s_mu_acv1(1:ns_mu_acv1));
        
        wQ0s_acv2(j) = mean(w(1:ns_Q0_acv2).*Q0s(1:ns_Q0_acv2));
        Q1s_acv2 = Q1(ss(1:ns_Q1_acv2))<0;
        wQ1s_acv2(j) = mean(w(1:ns_Q1_acv2).*Q1s_acv2(1:ns_Q1_acv2));
        ss_mu_acv2 = random(gm,ns_mu_acv2);
        Q1s_mu_acv2 = Q1(ss_mu_acv2)<0;
        w_mu_acv2 = mvnpdf(ss_mu_acv2,mu,std)./qce(ss_mu_acv2);
        mu_acv2(j) = (sum(w(1:ns_Q1_acv2).*Q1s_acv2(1:ns_Q1_acv2))+sum(w_mu_acv2(1:ns_mu_acv2).*Q1s_mu_acv2(1:ns_mu_acv2)))/(ns_Q1_acv2+ns_mu_acv2);
    end
    
    writematrix(wQ0s,append(path,'wQ0s.txt'));
    writematrix(wQ0s_cv,append(path,'wQ0s_cv.txt'));
    writematrix(wQ1s_cv,append(path,'wQ1s_cv.txt'));
    writematrix(wQ0s_acv1,append(path,'wQ0s_acv1.txt'));
    writematrix(wQ1s_acv1,append(path,'wQ1s_acv1.txt'));
    writematrix(mu_acv1,append(path,'mu_acv1.txt'));
    writematrix(wQ0s_acv2,append(path,'wQ0s_acv2.txt'));
    writematrix(wQ1s_acv2,append(path,'wQ1s_acv2.txt'));
    writematrix(mu_acv2,append(path,'mu_acv2.txt'));
    
%     m0 = mean(wQ0s);
%     m1 = mean(wQ1s);
%     m1_acv1 = mean(wQ1s_acv1);
%     m1_acv2 = mean(wQ1s_acv2);
%     me = m0+a*(m1-prob1);
%     me_acv1 = m0+a*(m1-m1_acv1);
%     me_acv2 = m0+a*(m1-m1_acv2);
%     
%     prob0
%     prob1
%     
%     display('prob0./m0')
%     prob0./m0
%     display('prob0./me')
%     prob0./me'
%     display('prob0./me_acv1')
%     prob0./me_acv1'
%     display('prob0./me_acv2')
%     prob0./me_acv2'

    v0 = var(wQ0s);
    cov_cv = cov(wQ0s_cv,wQ1s_cv);
    v0_cv = cov_cv(1,1);
    v1_cv = cov_cv(2,2);
    ve = v0_cv + a.^2*v1_cv + 2*a*cov_cv(1,2);
    
    cov_acv1 = cov(wQ0s_acv1,wQ1s_acv1);
    v0_acv1 = cov_acv1(1,1);
    v1_acv1 = cov_acv1(2,2);
    cov_Q1_mu_acv1 = cov(wQ1s_acv1,mu_acv1);
    v_mu_acv1 = cov_Q1_mu_acv1(2,2);
    ve_acv1 = v0_acv1 + a.^2*(v1_acv1+v_mu_acv1) + 2*a*cov_acv1(1,2);
    
    cov_acv2 = cov(wQ0s_acv2,wQ1s_acv2);
    v0_acv2 = cov_acv2(1,1);
    v1_acv2 = cov_acv2(2,2);
    cov_Q1_mu_acv2 = cov(wQ1s_acv2,mu_acv2);
    v_mu_acv2 = cov_Q1_mu_acv2(2,2);
    cov_0acv2 = cov(wQ0s_acv2,mu_acv2);
    ve_acv2 = v0_acv2 + a.^2*(v1_acv2+v_mu_acv2-2*cov_Q1_mu_acv2(1,2)) + 2*a*(cov_acv2(1,2)-cov_0acv2(1,2));
    
%     covar = cov(wQ0s,wQ1s);
%     covar_acv1 = cov(wQ1s,wQ1s_acv1);
%     covar_acv2 = cov(wQ1s,wQ1s_acv2);
%     covar_0acv2 = cov(wQ0s,wQ1s_acv2);
%     v0 = covar(1,1);
%     v1 = covar(2,2);
%     v_acv1 = covar_acv1(2,2);
%     v_acv2 = covar_acv2(2,2);
%     ve = v0 + a.^2*v1 + 2*a*covar(1,2);
%     ve_acv1 = v0 + a.^2*(v1+v_acv1) + 2*a*covar(1,2);
%     ve_acv2 = v0 + a.^2*(v1+v_acv2-2*covar_acv2(1,2)) + 2*a*(covar(1,2)-covar_0acv2(1,2));
      
    display('v./v0')
    v/v0
    display('v/ve')
    v./ve'
    display('ve./v0')
    ve'./v0
    display('ve_acv1./v0')
    ve_acv1'./v0
    display('ve_acv2./v0')
    ve_acv2'./v0
    
    min(ve'./v0)
    min(ve_acv1'./v0)
    min(ve_acv2'./v0)
    
%     a_opt = - covar(1,2)/v1
%     a_opt_acv1 = - covar(1,2)/(v1+v_acv1)
    
    set(0,'defaultLineLineWidth',0.7);
    set(0,'defaultLineMarkerSize',2);
    fontsize = 8;
    width = 2.5;
    height = 2.5;
    
    figs(1) = figure('Units','inches',...
        'Position',[0 0 width height],...
        'visible','off',...
        'PaperPositionMode','auto');
    hold on
    plot(a,ve,'-o',a,v0*ones(1,length(a)),'-o')
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times')
    legend({'$v_e$','$v_0$'},...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times',...
        'Location','NorthEast')
    ylabel('Variance',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times')
    xlabel({'$\alpha$'},...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times')
    fn(1) = "vev0.eps";
    hold off
    
    figs(2) = figure('Units','inches',...
        'Position',[0 0 width height],...
        'visible','on',...
        'PaperPositionMode','auto');
    hold on
%     plot(a,ve./v0,'-d',a,ve_acv1./v0,'-v',a,ve_acv2./v0,'-s',a,ones(1,length(a)),'-o')
    plot(a,ve./v0,'-d',a,ve_acv1./v0,'-v',a,ve_acv2./v0,'-s',a,ones(1,length(a)),'-o')
    yticks([1 3 5 7 9 10])
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times',...
        'Box','on')
    legend({'$v_e/v_0$','$v_1/v_0$','$v_2/v_0$'},...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times',...
        'Location','North')
    ylabel('Variance ratio',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times')
    xlabel({'$\alpha$'},...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times')
    fn(2) = "ve_v0_acv.eps";
    hold off
    
    q0 = @(x) ((Q0(x)<0).*mvnpdf(x,mu,std))/prob0;
    q1 = @(x) ((Q1(x)<0).*mvnpdf(x,mu,std))/prob1;
    x = linspace(2,5.5,1000)';
    figs(3) = figure('Units','inches',...
        'Position',[0 0 width height],...
        'visible','off',...
        'PaperPositionMode','auto');
    hold on
    plot(x,q0(x),'-',x,q1(x),'-',x,qce(x),'-')
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times','Box','on')
    legend({'$q_0(z)$','$q_1(z)$','$q_{CE}(z)$'},...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times',...
        'Location','NorthEast')
    xlabel({'$z$'},...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times')
    fn(3) = "qce.eps";
    hold off
    
    figs(4) = figure('Units','inches',...
        'Position',[0 0 width height],...
        'visible','off',...
        'PaperPositionMode','auto');
    hold on
    plot(a,log(ve),'-o',a,log(v0*ones(1,length(a))),'-o')
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times')
    legend({'$\ln(v_e)$','$\ln(v_0)$'},...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times',...
        'Location','NorthEast')
    ylabel('Logarithm of Variance',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times')
    xlabel({'$\alpha$'},...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times')
    fn(4) = "logvev0.eps";
    hold off
    
    for k = 1:length(figs)
        % print each figure in figs to a separate .eps file
        print(figs(k), '-depsc2', sprintf('%s', append(path,fn(k))))
    end
    
end