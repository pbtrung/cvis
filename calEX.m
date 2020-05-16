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
    writematrix(l*m,'C:\Users\Nathan\Data\U_lm.txt');
    writematrix(Umean,'C:\Users\Nathan\Data\Umean.txt');
    writematrix(Uvar,'C:\Users\Nathan\Data\Uvar.txt');
   
    m = max(Uy1(:));
    l = 0.05:0.05:1;
    U1mean(1:length(l)) = 0;
    U1var(1:length(l)) = 0;
    for j = 1:length(l)
        fprintf('iter: %d, l: %f, l*m: %f\n',j,l(j),l(j)*m);
        U1mean(j) = mean(mean(l(j)*m-Uy1<0,2));
        U1var(j) = var(mean(l(j)*m-Uy1<0,2));
    end
    writematrix(l'*m,'C:\Users\Nathan\Data\U1_lm_10x10_EX.txt');
    writematrix(U1mean','C:\Users\Nathan\Data\U1mean_10x10.txt');
    writematrix(U1var','C:\Users\Nathan\Data\U1var_10x10.txt');
    
end