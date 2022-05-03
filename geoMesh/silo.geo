h=5;
r1 = 30;
r2 = 40;

Point(1) = {0,0, 0, h};
Point(2) = {20,0, 0, h};
Point(3) = {60, 60, 0, h};
Point(4) = {60,260, 0, h};
Point(5) = {0,260, 0, h};

Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,5};
Line(5)  = {5,1};
Line Loop(5) = {1, 2, 3, 4,5};

Plane Surface(1) = {5};
