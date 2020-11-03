function benchmark()
    
    format long;
    rng('default');
    
    [filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
    repath = fullfile(filepath,'results');
    
    t0 = readmatrix(fullfile(repath,'t_30x30_midNode.txt'));
    t1 = readmatrix(fullfile(repath,'t_10x10_midNode.txt'));
    
    disp(mean(t0)/mean(t1));
    
end