function calEX(p,p1)

    Uy = readmatrix(p);
    Uy1 = readmatrix(p1);
    
    m = max(Uy(:));
    l = 0.05:0.05:1;
    Umean(1:length(l)) = 0;
    Uvar(1:length(l)) = 0;
    for j = 1:length(l)
        fprintf('iter: %d, l: %f, l*m: %f\n',j,l(j),l(j)*m);
        Umean(j) = mean(mean(l(j)*m-Uy<0,2));
        Uvar(j) = var(mean(l(j)*m-Uy<0,2));
    end
    writematrix(l*m,'U_lm.txt');
    writematrix(Umean,'Umean.txt');
    writematrix(Uvar,'Uvar.txt');
   
    m = max(Uy1(:));
    l = 0.05:0.05:1;
    U1mean(1:length(l)) = 0;
    U1var(1:length(l)) = 0;
    for j = 1:length(l)
        fprintf('iter: %d, l: %f, l*m: %f\n',j,l(j),l(j)*m);
        U1mean(j) = mean(mean(l(j)*m-Uy1<0,2));
        U1var(j) = var(mean(l(j)*m-Uy1<0,2));
    end
    writematrix(l*m,'U1_lm.txt');
    writematrix(U1mean,'U1mean.txt');
    writematrix(U1var,'U1var.txt');
    
end

