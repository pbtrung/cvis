function [cov_Q0Q1, rho] = cvis_ce_rf_read_cov()
    
    format long;
    
    [filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
    repath = fullfile(filepath,'Example 2');

    wQ0s_MC = readmatrix(fullfile(repath,'wQ0s_MC_cov.txt'));
    wQ1s_MC = readmatrix(fullfile(repath,'wQ1s_MC_cov.txt'));
    
    cov_Q0Q1_MC = cov(wQ0s_MC,wQ1s_MC);
    cov_Q0Q1 = cov_Q0Q1_MC(1,2)
    rho = cov_Q0Q1_MC(1,2)/sqrt(cov_Q0Q1_MC(1,1)*cov_Q0Q1_MC(2,2))
    
end