function cvis_rs(path)
    
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
    n_MC = 100000;
    
    tol = 1.2;
    prob2 = 1-normcdf(l1-tol);
    Q2 = @(x) l1-tol-x;
    q = @(x) ((Q2(x)<0).*normpdf(x))/prob2;
    
    % rejection sampling
    s_MC = rs(q,n_MC);
    Q0s_MC = Q0(s_MC(:))<0;
    w_MC = normpdf(s_MC)./q(s_MC);
    wQ0s_MC = w_MC.*Q0s_MC;
    
    % c1 = c*c2
    c = 10;
    n0 = fix(n_MC*c/(c+1));
    n1 = n0;
    
    Q0s = Q0s_MC(1:n0);
    w0 = normpdf(s_MC(1:n0))./q(s_MC(1:n0));
    wQ0s = w0.*Q0s;
    Q1s = Q1(s_MC(1:n1))<0;
    w1 = normpdf(s_MC(1:n1))./q(s_MC(1:n1));
    wQ1s = w1.*Q1s;
    
    m_MC = mean(wQ0s_MC);
    m0 = mean(wQ0s);
    m1 = mean(wQ1s);
    m_CV = m0 + a*(m1-prob1);
    
    disp('prob0/m_MC');
    disp(prob0/m_MC);
    disp('prob0/m0');
    disp(prob0/m0);
    disp('prob0/m_CV');
    disp(prob0./m_CV');
    
    v0_MC = var(wQ0s_MC)/n_MC;
    cov_Q0Q1 = cov(wQ0s,wQ1s);
    cov_MC_Q0Q1 = cov_Q0Q1/n0;
    v_MC_Q0 = cov_MC_Q0Q1(1,1);
    v_MC_Q1 = cov_MC_Q0Q1(2,2);
    v_CV = v_MC_Q0 + a.^2*v_MC_Q1 + 2*a*cov_MC_Q0Q1(1,2);
    
    disp('v_CV/v0_MC');
    disp(v_CV'/v0_MC);
    
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
    plot(a,v_CV'/v0_MC,'-s',a,ones(1,length(a)),'-o')
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times','Box','on')
    legend({'$Var_{\hat{q}}[\mathcal{Q}^{MF}_{\hat{q},n}]$','$Var_{\hat{q}}[\mathcal{Q}^{MFIS}_{\hat{q},n}]$'},...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times',...
        'Location','Best')
    ylim([0.45 2.75]);
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
    fn(1) = "v_CV__v0_MC.eps";
    hold off
    
    for k = 1:length(figs)
        % print each figure in figs to a separate .eps file
        print(figs(k), '-depsc2', sprintf('%s', append(path,fn(k))))
    end
    
end

function s = rs(q,n)
    umin = 0;
    umax = 7;
    rnsamples = n*50;
    u = umin + (umax-umin)*rand(rnsamples,1);
    sample_value = q(u);
    max_value = max(sample_value);
    accepted = rand(rnsamples,1) < (sample_value/max_value);
    s = u(accepted,:);
    if length(s) < n
        error('length(s) < n');
    else
        s = s(1:n);
    end
end