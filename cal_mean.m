function cal_mean(p0,p1)

    Uy0 = readmatrix(sprintf('%s',p0));
%     m = max(Uy0(:));
%     l = 0.97*m;
    l = 118.923;
    mean(l-Uy0<0) % 0.001173
    var(l-Uy0<0)
    
    Uy1 = readmatrix(sprintf('%s',p1));
%     m = max(Uy1(:));
%     l = 0.9*m
    l = 108.510;
    mean(l-Uy1<0) % 0.022428
    var(l-Uy1<0)
    
end