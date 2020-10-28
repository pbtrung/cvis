function testKE()
    KE = ElmStiffnessMatrix(60,20,60,20);
    KE9 = KE99();
    
    if isequal(ismembertol(KE,KE9),true(8,8)) == false
        disp('NOT EQUAL');
    else
        disp('EQUAL');
    end
    
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