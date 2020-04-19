function [Psi,lambda] = KL(lx,ly,Lx,Ly,nelx,nely)

neig = 10;

K1 = @(s,t) exp(-abs(s-t)/Lx);
F1 = chebop(@(u) fred(K1, u));
[Psi1,Lambda1] = eigs(F1,neig,'lm');
Psi1 = chebfun(Psi1);
[lambda1,idx1] = sort(diag(Lambda1),'descend');
Psi1 = Psi1(:,idx1);

K2 = @(s,t) exp(-abs(s-t)/Ly);
F2 = chebop(@(u) fred(K2, u));
[Psi2,Lambda2] = eigs(F2,neig,'lm');
Psi2 = chebfun(Psi2);
[lambda2,idx2] = sort(diag(Lambda2),'descend');
Psi2 = Psi2(:,idx2);

lelx = lx/nelx;
lely = ly/nely;
x1 = 0:lelx/2:lx;
x1 = x1(2:2:end-1);
x2 = 0:lely/2:ly;
x2 = x2(2:2:end-1);

lbd(1:neig,1:neig)= 0;
for i = 1:neig
    lbd(i,:) = lambda1(i)*lambda2;
end
[lam_idx1,lam_idx2] = ndgrid(1:size(lbd,1),1:size(lbd,2));
[lambda,idx] = sort(lbd(:),'descend');

neig = neig*2;
lambda = lambda(1:neig);
% L = diag(sqrt(lambda(1:neig)));
lam_idx1 = lam_idx1(idx);
lam_idx2 = lam_idx2(idx);
Psi1 = Psi1(:,lam_idx1(1:neig));
Psi2 = Psi2(:,lam_idx2(1:neig));

Psi(1:nely,1:nelx,1:neig) = 0;
for i = 1:nelx
    for j = 1:nely
        Psi(j,i,:) = Psi1(x1(i)).*Psi2(x2(j));
    end
end

% for i = 1:nelx
%     for j = 1:nely
%         Z = randn(neig,1);
%         X = squeeze(Psi(j,i,:))'*(L*Z);
%         E(j,i) = a+(b-a)*normcdf(X);
%     end
% end
% E

end