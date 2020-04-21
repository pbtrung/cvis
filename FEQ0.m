% FE-ANALYSIS
function Udof = FEQ0(lx,ly,nelx,nely,dof,force,x,P,L)
    
    a = 1;
    b = 2;
    LZ = L*x;
    PZ = dot(P',LZ);
    X = reshape(PZ,nelx,nely)';
    E = a+(b-a)*normcdf(X);

    K = sparse(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
    F = sparse(2*(nely+1)*(nelx+1),1);
    U = zeros(2*(nely+1)*(nelx+1),1);
    KE = ElmStiffnessMatrix(lx,ly,nelx,nely);
    
    for elx = 1:nelx
        for ely = 1:nely
            n1 = (nely+1)*(elx-1) + ely;
            n2 = (nely+1)*elx + ely;
            edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
            K(edof,edof) = K(edof,edof) + E(ely,elx)*KE;
        end
    end
    % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
    F(dof,1) = force;
    fixeddofs   = 1:2*(nely+1);
    alldofs     = 1:2*(nely+1)*(nelx+1);
    freedofs    = setdiff(alldofs,fixeddofs);
    % SOLVING
    U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
    U(fixeddofs,:)= 0;
    Udof = U(dof);
end

function KE = ElmStiffnessMatrix(lx,ly,nelx,nely)
    E = 1.;
    nu = 0.3;
    h = 1;
    a = lx/nelx;
    b = ly/nely;
    c = b/a;
    
    ka = (2/c)*(1-nu);
    kb = (3/2)*(1+nu);
    kc = (3/2)*(1-3*nu);
    kd = (2*c)*(1-nu);
    KE = (E*h)/(12*(1-nu^2))*[ ...
        4*c+ka kb -4*c+ka/2 -kc -2*c-ka/2 -kb 2*c-ka kc
        kb 4/c+kd kc 2/c-kd -kb -2/c-kd/2 -kc -4/c+kd/2
        -4*c+ka/2 kc 4*c+ka -kb 2*c-ka -kc -2*c-ka/2 kb
        -kc 2/c-kd -kb 4/c+kd kc -4/c+kd/2 kb -2/c-kd/2
        -2*c-ka/2 -kb 2*c-ka kc 4*c+ka kb -4*c+ka/2 -kc
        -kb -2/c-kd/2 -kc -4/c+kd/2 kb 4/c+kd kc 2/c-kd
        2*c-ka -kc -2*c-ka/2 kb -4*c+ka/2 kc 4*c+ka -kb
        kc -4/c+kd/2 kb -2/c-kd/2 -kc 2/c-kd -kb 4/c+kd
        ];
end