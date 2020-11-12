function cvis_ce_acv_out()

    format long;
    rng('default');
                    
    [filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
    repath = fullfile(filepath,'results','Example 3');
    
    wQ0s_MC = readmatrix(fullfile(repath,'wQ0s_MC.txt'));
    wQ1s_MC = readmatrix(fullfile(repath,'wQ1s_MC.txt'));
        
    nbatches = 1000;
    
    % c0 = c*c1
    c = 37;
    
    %%
    n_MC = 4*10^5;
    v0_MC = var(wQ0s_MC(1:n_MC))/n_MC;
    
    n0_CV = fix(n_MC*c/(c+1));
    n1_CV = n0_CV;
    
    n0_CV_EN = fix(n0_CV/nbatches);
    n1_CV_EN = fix(n1_CV/nbatches);
    wQ0s_MC = wQ0s_MC(1:n0_CV_EN*nbatches);
    wQ1s_MC = wQ1s_MC(1:n1_CV_EN*nbatches);
    wQ0s_MC = reshape(wQ0s_MC,n0_CV_EN,nbatches);
    wQ1s_MC = reshape(wQ1s_MC,n1_CV_EN,nbatches);
    
    wQ0s_CV_EN(1:nbatches) = 0;
    wQ1s_CV_EN(1:nbatches) = 0;
    for i = 1:nbatches
        wQ0s_CV_EN(i) = mean(wQ0s_MC(:,i));
        wQ1s_CV_EN(i) = mean(wQ1s_MC(:,i));
    end
    
    cov_Q0Q1_CV_EN = cov(wQ0s_CV_EN,wQ1s_CV_EN);
    cov_MC_Q0Q1_CV_EN = cov_Q0Q1_CV_EN/nbatches;
    v_MC_Q0_CV_EN = cov_MC_Q0Q1_CV_EN(1,1);
    v_MC_Q1_CV_EN = cov_MC_Q0Q1_CV_EN(2,2);
    astar = -cov_MC_Q0Q1_CV_EN(1,2)/v_MC_Q1_CV_EN;
    v_CV_EN = v_MC_Q0_CV_EN + astar^2*v_MC_Q1_CV_EN + 2*astar*cov_MC_Q0Q1_CV_EN(1,2);
    disp('v_CV_EN/v0_MC');
    disp(v_CV_EN/v0_MC);
    
    %%
    mu1_IS = readmatrix(fullfile(repath,'wQ1s_IS.txt'));
    
    r = 3.5;
    n0_IS = fix(n_MC*c/(c+r+1));
    n1_IS = n0_IS;
    m1_IS = fix(r*n0_IS);
    
    n0_IS_EN = fix(n0_IS/nbatches);
    n1_IS_EN = fix(n1_IS/nbatches);
    m1_IS_EN = fix(m1_IS/nbatches);
    
    wQ0s_IS = wQ0s_MC(1:n0_IS_EN*nbatches);
    wQ1s_IS = wQ1s_MC(1:n1_IS_EN*nbatches);
    mu1_IS = mu1_IS(1:m1_IS_EN*nbatches);
    wQ0s_IS = reshape(wQ0s_IS,n0_IS_EN,nbatches);
    wQ1s_IS = reshape(wQ1s_IS,n1_IS_EN,nbatches);
    mu1_IS = reshape(mu1_IS,m1_IS_EN,nbatches);
    
    wQ0s_IS_EN(1:nbatches) = 0;
    wQ1s_IS_EN(1:nbatches) = 0;
    mu1_IS_EN(1:nbatches) = 0;
    for i = 1:nbatches
        wQ0s_IS_EN(i) = mean(wQ0s_IS(:,i));
        wQ1s_IS_EN(i) = mean(wQ1s_IS(:,i));
        mu1_IS_EN(i) = (sum(wQ1s_IS(:,i)) + sum(mu1_IS(:,i)))/(n1_IS_EN + m1_IS_EN);
    end
    
    cov_Q0Q1_IS_EN = cov(wQ0s_IS_EN,wQ1s_IS_EN);
    cov_MC_Q0Q1_IS_EN = cov_Q0Q1_IS_EN/nbatches;
    v_MC_Q0_IS_EN = cov_MC_Q0Q1_IS_EN(1,1);
    v_MC_Q1_IS_EN = cov_MC_Q0Q1_IS_EN(2,2);
    cov_MC_Q0mu1_IS_EN = cov(wQ0s_IS_EN,mu1_IS_EN)/nbatches;
    v_MC_mu1_IS_EN = cov_MC_Q0mu1_IS_EN(2,2);
    cov_MC_Q1mu1_IS_EN = cov(wQ1s_IS_EN,mu1_IS_EN)/nbatches;
    astar = -(cov_MC_Q0Q1_IS_EN(1,2)-cov_MC_Q0mu1_IS_EN(1,2))...
        /(v_MC_Q1_IS_EN + v_MC_mu1_IS_EN - 2*cov_MC_Q1mu1_IS_EN(1,2));
    v_IS_EN = v_MC_Q0_IS_EN +...
        astar^2*(v_MC_Q1_IS_EN + v_MC_mu1_IS_EN - 2*cov_MC_Q1mu1_IS_EN(1,2))...
        + 2*astar*(cov_MC_Q0Q1_IS_EN(1,2)-cov_MC_Q0mu1_IS_EN(1,2));
    disp('v_IS_EN/v0_MC');
    disp(v_IS_EN/v0_MC);

end