function cvis_ce_rf_acv()
    
    format long;
    rng('default');
                
    a = linspace(-1.5,0.5,33);
    
    [filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
    repath = fullfile(filepath,'results');
    
    n_MC = 4*10^5;
    wQ0s_MC = readmatrix(fullfile(repath,'wQ0s_MC.txt'));
    wQ1s_MC = readmatrix(fullfile(repath,'wQ1s_MC.txt'));
    
    % c0 = c*c1
    c = 11;
    
    %%
    n0_CV = fix(n_MC*c/(c+1));
    n1_CV = n0_CV;
    
    wQ0s_CV = wQ0s_MC(1:n0_CV);
    wQ1s_CV = wQ1s_MC(1:n1_CV);
    
    v0_MC = var(wQ0s_MC(1:n_MC))/n_MC;
    cov_Q0Q1_CV = cov(wQ0s_CV,wQ1s_CV);
    cov_MC_Q0Q1_CV = cov_Q0Q1_CV/n0_CV;
    v_MC_Q0_CV = cov_MC_Q0Q1_CV(1,1);
    v_MC_Q1_CV = cov_MC_Q0Q1_CV(2,2);
    v_CV = v_MC_Q0_CV + a.^2*v_MC_Q1_CV + 2*a*cov_MC_Q0Q1_CV(1,2);
    
    %%
    r = 3;
    n0_IS = fix(n_MC*c/(c+r+1));
    n1_IS = n0_IS;
    m1_IS = fix(r*n0_IS);
    r1 = (n1_IS+m1_IS)/n0_IS;
    
    wQ0s_IS = wQ0s_MC(1:n0_IS);
    wQ1s_IS = wQ1s_MC(1:n1_IS);
        
    %%
    cov_Q0Q1_IS = cov(wQ0s_IS,wQ1s_IS);
    cov_MC_Q0Q1_IS = cov_Q0Q1_IS/n0_IS;
    v_MC_Q0_IS = cov_MC_Q0Q1_IS(1,1);
    v_MC_Q1_IS = cov_MC_Q0Q1_IS(2,2);
    v_IS = v_MC_Q0_IS + a.^2*(v_MC_Q1_IS*((r1-1)/r1)) + 2*a*(cov_MC_Q0Q1_IS(1,2)*((r1-1)/r1));
    
    disp(min(v_CV'/v0_MC));
    disp(min(v_IS'/v0_MC));
    
    %%
    set(0,'defaultLineLineWidth',0.7);
    set(0,'defaultLineMarkerSize',2);
    fontsize = 8;
    width = 2.5;
    height = 2.5;
    
    figs(1) = figure('Units','inches',...
        'Position',[0 0 width height],...
        'visible','on',...
        'PaperPositionMode','auto');
    hold on
    plot(a,v_CV'/v0_MC,'-d',a,v_IS'/v0_MC,'-v',a,ones(1,length(a)))
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times','Box','on')
    legend({'$v_{CV}$','$v_{IS}$','$v_0$'},...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times',...
        'Location','Best')
    ylim([0.35 3.1]);
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
    fn(1) = "v_ACV__v0_MC.eps";
    hold off
    
    for k = 1:length(figs)
        % print each figure in figs to a separate .eps file
        print(figs(k), '-depsc2', fullfile(repath,fn(k)))
    end
    
end