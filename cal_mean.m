function cal_mean()

    [filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
    repath = fullfile(filepath,'results');

    Uy0 = readmatrix(fullfile(repath,'Uy0.txt'));
%     m = max(Uy0(:));
%     l = 0.97*m;
    l = 118.923;
    disp(mean(l-Uy0<0)) 
    % 0.001173
    disp(var(l-Uy0<0)/length(Uy0))
    % 1.171625242625240e-09
    
    Uy1 = readmatrix(fullfile(repath,'Uy1.txt'));
%     m = max(Uy1(:));
%     l = 0.9*m
    l = 108.510;
    disp(mean(l-Uy1<0)) 
    % 0.022428
    disp(var(l-Uy1<0)/length(Uy1))
    % 2.192500674100718e-08
    
end