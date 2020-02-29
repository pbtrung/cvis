% FE-ANALYSIS
function U = FE3(nelx,nely,dof,KE,force)
    nf = length(force);
    K = sparse(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
    F = sparse(2*(nely+1)*(nelx+1),nf); 
    U = zeros(2*(nely+1)*(nelx+1),nf);
    for elx = 1:nelx
      for ely = 1:nely
        n1 = (nely+1)*(elx-1) + ely; 
        n2 = (nely+1)*elx + ely;
        edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
        K(edof,edof) = K(edof,edof) + KE;
      end
    end
    % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
    fixeddofs   = 1:2*(nely+1);
    alldofs     = 1:2*(nely+1)*(nelx+1);
    freedofs    = setdiff(alldofs,fixeddofs);
    % SOLVING
    F(dof,:) = force;
    for j = 1:nf
        U(freedofs,j) = K(freedofs,freedofs) \ F(freedofs,j);
        U(fixeddofs,j)= 0;
    end
end