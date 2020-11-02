E = 10000;
poisson = 0.30;
kapa = 5/6;
L = 1;
nelx = 20;
nely = 20;
h = 20*0.1;
x = [1 1 1 1 h h h h];
FE_plate(nelx,nely,E,poisson,kapa,L,x)