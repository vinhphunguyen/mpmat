function nodes=getNodesForParticle2D(x, y, dx, dy, numx)

xi = floor ( x/dx ) + 1;
yi = floor ( y/dy )    ;

n1 = xi + (numx+1)*yi;
n4 = xi + (numx+1)*(yi+1);

n2 = n1 + 1;
n3 = n4 + 1;

nodes(1)   = n1;
nodes(2) = n2;
nodes(3) = n3;
nodes(4) = n4;
