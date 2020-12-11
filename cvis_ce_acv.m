function cvis_ce_acv()

    format long;
    rng('default');

    mu = 0;
    std = 1;
    w0 = 3;
    w1 = 2.8;
    l0 = mu + w0*std;
    l1 = mu + w1*std;
    
    Q0 = @(x) l0-x;
    Q1 = @(x) l1-x;
    
    prob0 = 1-normcdf(l0);
    prob1 = 1-normcdf(l1);
    
    a = linspace(-1.5,0.5,33);
    
    % definition of the random variables
    d      = 1;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);
    
    % CE method
    N      = 3000;    % total number of samples for each level
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
    
    % limit state function
    g = @(x) Q1(x);
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    qce = @(x) pdf(gm,x);
    
    n_MC = 5*10^5;
    s_MC = random(gm, n_MC);
    Q0s_MC = Q0(s_MC(:))<0;
    w_MC = mvnpdf(s_MC,mu,std)./qce(s_MC);
    wQ0s_MC = w_MC.*Q0s_MC;
    
    % c1 = c*c2
    c = 30;
    
    %%
    n0_CV = fix(n_MC*c/(c+1));
    n1_CV = n0_CV;
    
    Q0s_CV = Q0s_MC(1:n0_CV);
    w0_CV = mvnpdf(s_MC(1:n0_CV),mu,std)./qce(s_MC(1:n0_CV));
    wQ0s_CV = w0_CV.*Q0s_CV;
    Q1s_CV = Q1(s_MC(1:n1_CV))<0;
    w1_CV = mvnpdf(s_MC(1:n1_CV),mu,std)./qce(s_MC(1:n1_CV));
    wQ1s_CV = w1_CV.*Q1s_CV;
    
    v0_MC = var(wQ0s_MC)/n_MC;
    cov_Q0Q1_CV = cov(wQ0s_CV,wQ1s_CV);
    cov_MC_Q0Q1_CV = cov_Q0Q1_CV/n0_CV;
    v_MC_Q0_CV = cov_MC_Q0Q1_CV(1,1);
    v_MC_Q1_CV = cov_MC_Q0Q1_CV(2,2);
    astar = -cov_MC_Q0Q1_CV(1,2)/v_MC_Q1_CV;
    v_CV = v_MC_Q0_CV + a.^2*v_MC_Q1_CV + 2*a*cov_MC_Q0Q1_CV(1,2);
    v_CV_min = v_MC_Q0_CV + astar^2*v_MC_Q1_CV + 2*astar*cov_MC_Q0Q1_CV(1,2);
    disp('v_CV_min/v0_MC');
    disp(v_CV_min/v0_MC);
    
%     m0_CV = mean(wQ0s_CV);
%     m1_CV = mean(wQ1s_CV);
%     m_CV = m0_CV + a*(m1_CV-prob1);
%     disp('prob0/m_CV');
%     disp(abs(prob0-m_CV')/prob0*100);
%     % from cvis_truth
%     m_naive = 0.00134859;
%     disp('prob0/m_naive');
%     disp(abs(prob0-m_naive)/prob0*100);
    
    %%
    r = 3.5;
    n0_IS = fix(n_MC*c/(c+r+1));
    n1_IS = n0_IS;
    m1_IS = fix(r*n0_IS);
    r1 = (n1_IS+m1_IS)/n0_IS;
    
    Q0s_IS = Q0s_MC(1:n0_IS);
    w0_IS = mvnpdf(s_MC(1:n0_IS),mu,std)./qce(s_MC(1:n0_IS));
    wQ0s_IS = w0_IS.*Q0s_IS;
    Q1s_IS = Q1(s_MC(1:n1_IS))<0;
    w1_IS = mvnpdf(s_MC(1:n1_IS),mu,std)./qce(s_MC(1:n1_IS));
    wQ1s_IS = w1_IS.*Q1s_IS;
    
    %%
    cov_Q0Q1_IS = cov(wQ0s_IS,wQ1s_IS);
    cov_MC_Q0Q1_IS = cov_Q0Q1_IS/n0_IS;
    v_MC_Q0_IS = cov_MC_Q0Q1_IS(1,1);
    v_MC_Q1_IS = cov_MC_Q0Q1_IS(2,2);
    v_IS = v_MC_Q0_IS + a.^2*(v_MC_Q1_IS*((r1-1)/r1)) + 2*a*(cov_MC_Q0Q1_IS(1,2)*((r1-1)/r1));
    
%     m0s_IS = mean(wQ0s_IS);
%     m1s_IS = mean(wQ1s_IS);
%     mu1s_IS = (sum(wQ1s_IS) + sum(wmu_IS))/(n1_IS + m1_IS);
%     m_IS = m0s_IS + a*(m1s_IS-mu1s_IS);
%     disp('prob0/m_IS');
%     disp(abs(prob0-m_IS')/prob0*100);
    
%     disp('v_CV/v0_MC');
%     disp(v_CV'/v0_MC);
%     disp('v_IS/v0_MC');
%     disp(v_IS'/v0_MC);
    
    disp("min(v_CV'/v0_MC)");
    disp(min(v_CV'/v0_MC));
    disp("min(v_IS'/v0_MC)");
    disp(min(v_IS'/v0_MC));
    % from cvis_truth
    v_naive = 1.346771318479568e-11;
    disp(v_naive/v0_MC);
    disp(max(v_naive./v_CV'));
    disp(max(v_naive./v_IS'));
        
    %%
    set(0,'defaultLineLineWidth',0.7);
    set(0,'defaultLineMarkerSize',2);
    fontsize = 8;
    width = 2.5;
    height = 2.5;
    
    % from cvis_ce_acv_en nbatches = 1000
    [v_CV_EN__MC,v_IS_EN__MC] = cvis_ce_acv_en_out(c,r,1000);
    v_CV_EN__MC_3for = 0.763171802502491;
    v_IS_EN__MC_3for = 0.934273161041501;
    
    figs(1) = figure('Units','inches',...
        'Position',[0 0 width height],...
        'visible','on',...
        'PaperPositionMode','auto');
    hold on
    plot(a,v_CV'/v0_MC,'-d',a,v_IS'/v0_MC,'-v',a,ones(1,length(a)),...
        a,v_CV_EN__MC*ones(1,length(a)),a,v_IS_EN__MC*ones(1,length(a)),...
        a,v_CV_EN__MC_3for*ones(1,length(a)),a,v_IS_EN__MC_3for*ones(1,length(a)))
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times','Box','on')
    legend({'$v_{CV}$','$v_{IS}$','$v_0$','$\bar{v}_{CV}$','$\bar{v}_{IS}$',...
        '$\tilde{v}_{CV}$','$\tilde{v}_{IS}$'},...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times',...
        'Location','Best')
    ylim([0.1 2.5]);
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
    fn(1) = "v_ACV__v0_MC_30.eps";
    hold off
    
%     figs(2) = figure('Units','inches',...
%         'Position',[0 0 width height],...
%         'visible','on',...
%         'PaperPositionMode','auto');
%     hold on
%     plot(a,abs(prob0-m_CV')/prob0*100,'-d',a,abs(prob0-m_IS')/prob0*100,'-v',a,abs(prob0-m_naive)/prob0*100*ones(1,length(a)),'-o')
%     set(gca,...
%         'Units','normalized',...
%         'FontUnits','points',...
%         'FontWeight','normal',...
%         'FontSize',fontsize,...
%         'FontName','Times','Box','on')
%     legend({'$v_{CV}$','$v_{IS}$','$v_0$'},...
%         'interpreter','latex',...
%         'FontSize',fontsize,...
%         'FontName','Times',...
%         'Location','Best')
%     ylabel('Variance ratio',...
%         'FontUnits','points',...
%         'interpreter','latex',...
%         'FontSize',fontsize,...
%         'FontName','Times')
%     xlabel({'$\alpha$'},...
%         'FontUnits','points',...
%         'interpreter','latex',...
%         'FontWeight','normal',...
%         'FontSize',fontsize,...
%         'FontName','Times')
%     fn(2) = "v_CV__v0_MC.eps";
%     hold off

    q0 = @(x) ((Q0(x)<0).*mvnpdf(x,mu,std))/prob0;
    q1 = @(x) ((Q1(x)<0).*mvnpdf(x,mu,std))/prob1;
    x = linspace(2,5.5,1000)';
    figs(2) = figure('Units','inches',...
        'Position',[0 0 width height],...
        'visible','on',...
        'PaperPositionMode','auto');
    hold on
    plot(x,q0(x),'-',x,q1(x),'-',x,qce(x),'-')
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times','Box','on')
    legend({'$q_0(z)$','$q_1(z)$','$\hat{q}(z)$'},...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times',...
        'Location','Best')
    ylabel('',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times')
    xlabel({'$z$'},...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times')
    fn(2) = 'qce.eps';
    hold off
    
    [filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
    repath = fullfile(filepath,'results');
    
    for k = 1:length(figs)
        % print each figure in figs to a separate .eps file
        print(figs(k), '-depsc2', fullfile(repath,fn(k)))
    end
    
end