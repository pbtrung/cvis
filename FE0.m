function Udof = FE0(lx,ly,nelx,nely,dof,force,x,P,L,lower,upper)
    
    nx = length(x); % number of columns
    Udof(1:nx,1) = 0;

    KE = ElmStiffnessMatrix(lx,ly,nelx,nely);
    
    % DEFINE LOADS AND SUPPORTS
    fixeddofs = 1:2*(nely+1);
    alldofs   = 1:2*(nely+1)*(nelx+1);
    freedofs  = setdiff(alldofs,fixeddofs);
        
    parfor j = 1:nx
        LZ = L*x(:,j);
        PZ = dot(P',repmat(LZ,1,nelx*nely));
        X = reshape(PZ,nelx,nely)';
        E = lower + (upper-lower)*normcdf(X);
        
        F = sparse(2*(nely+1)*(nelx+1),1);
        F(dof,1) = force;
        K = sparse(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
        U = zeros(2*(nely+1)*(nelx+1),1);
        for elx = 1:nelx
            for ely = 1:nely
                n1 = (nely+1)*(elx-1) + ely;
                n2 = (nely+1)*elx + ely;
                edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
                K(edof,edof) = K(edof,edof) + E(ely,elx)*KE;
            end
        end
        
        % SOLVING
        U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
        U(fixeddofs,:)= 0;
        Udof(j) = U(dof);
    end
    
end