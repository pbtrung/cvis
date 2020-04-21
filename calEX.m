function calEX(p,p1)
    Uy = readmatrix(p);
    Uy1 = readmatrix(p1);
    
    m = max(Uy(:));
%     l = 0.99*m
    l = 84.104;
    mean(mean(l-Uy<0,2))
    var(mean(l-Uy<0,2))
    
    m = max(Uy1(:));
%     l = 0.98*m
    l = 118.351;
    mean(mean(l-Uy1<0,2))
    var(mean(l-Uy1<0,2))
end

