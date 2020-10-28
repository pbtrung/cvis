function U = FE(lx, ly, nelx, nely, dof, force, E, fl)
    K = sparse(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
    F = sparse(2*(nely+1)*(nelx+1),1);
    U = zeros(2*(nely+1)*(nelx+1),1);
    KE = ElmStiffnessMatrix(lx,ly,nelx,nely);
    
    for elx = 1:nelx
        for ely = 1:nely
            n1 = (nely+1)*(elx-1) + ely;
            n2 = (nely+1)*elx + ely;
            edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
            if fl == 0
                K(edof,edof) = K(edof,edof) + E(ely,elx)*KE;
            elseif fl == 1
                K(edof,edof) = K(edof,edof) + E*KE;
            end
        end
    end
    % DEFINE LOADS AND SUPPORTS
    F(dof,1)    = force;
    fixeddofs   = 1:2*(nely+1);
    alldofs     = 1:2*(nely+1)*(nelx+1);
    freedofs    = setdiff(alldofs,fixeddofs);
    % SOLVING
    U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
    U(fixeddofs,:)= 0;
end