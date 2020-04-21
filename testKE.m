function testKE()
    KE = ElmStiffnessMatrix(60,20,60,20);
    KE9 = KE99();
    
    if isequal(ismembertol(KE,KE9),true(8,8)) == false
        display('NOT EQUAL')
    else
        display('EQUAL')
    end
    
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

function KE = KE99()
    E = 1.; 
    nu = 0.3;
    k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
       -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
    KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                      k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                      k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                      k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                      k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                      k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                      k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                      k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end