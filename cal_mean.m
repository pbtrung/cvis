function cal_mean(p)

    Uy = readmatrix(sprintf('%s',p));
    
    m = max(Uy(:));
%     l = 0.97*m
    l = 118.923;
    mean(l-Uy<0)
    var(l-Uy<0)

end