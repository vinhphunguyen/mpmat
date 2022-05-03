h=10;

Point(1)={0,0,0};
Point(2)={20000,0,0};
Point(3)={20000,10000,0};
Point(4)={10000,10000,0};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};

Line Loop(1)={1,2,3,4};

Plane Surface(1)={1};

Physical Surface(100)={1};
Physical Line(200)={1};
Physical Line(210)={2};


