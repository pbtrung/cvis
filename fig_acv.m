function fig_acv(save_path,data_path,p_wQ0s,p_wQ0s_cv,p_wQ1s_cv,p_wQ0s_acv1,p_wQ1s_acv1,p_mu_acv1,p_wQ0s_acv2,p_wQ1s_acv2,p_mu_acv2)

    wQ0s = readmatrix(append(data_path,p_wQ0s));
    
    wQ0s_cv = readmatrix(append(data_path,p_wQ0s_cv));
    wQ1s_cv = readmatrix(append(data_path,p_wQ1s_cv));
    
    wQ0s_acv1 = readmatrix(append(data_path,p_wQ0s_acv1));
    wQ1s_acv1 = readmatrix(append(data_path,p_wQ1s_acv1));
    mu_acv1 = readmatrix(append(data_path,p_mu_acv1));
    
    wQ0s_acv2 = readmatrix(append(data_path,p_wQ0s_acv2));
    wQ1s_acv2 = readmatrix(append(data_path,p_wQ1s_acv2));
    mu_acv2 = readmatrix(append(data_path,p_mu_acv2));
    
    a = linspace(-2.3,0.2,40);
    
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
        
    set(0,'defaultLineLineWidth',0.7);
    set(0,'defaultLineMarkerSize',2);
    fontsize = 8;
    width = 2.5;
    height = 2.5;
    
    ve'./v0
    ve_acv1'./v0
    ve_acv2'./v0
    
    min(ve'./v0)
    min(ve_acv1'./v0)
    min(ve_acv2'./v0)
    
    figs(1) = figure('Units','inches',...
        'Position',[0 0 width height],...
        'visible','on',...
        'PaperPositionMode','auto');
    hold on
    plot(a,ve./v0,'-d',a,ve_acv1./v0,'-v',a,ve_acv2./v0,'-s',a,ones(1,length(a)),'-o');
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times','Box','on')
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
    fn(1) = "ve_v0_acv.eps";
    hold off
    
    for k = 1:length(figs)
        % print each figure in figs to a separate .eps file  
        print(figs(k), '-depsc2', sprintf('%s', append(save_path,fn(k))))
    end
    
end