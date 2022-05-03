h=40;
h1=h/6;
h10=h1/3;
h2=h/4;


l=200;
a=5;
b=25;

Point(1)={0,0,0,h};
Point(2)={l,0,0,h};
Point(3)={l,l,0,h};
Point(4)={0,l,0,h};
Point(5)={l,0.5*l-0.5*a,0,h1};
Point(6)={l-b,0.5*l-0.5*a,0,h10};
Point(7)={l-b,0.5*l+0.5*a,0,h10};
Point(8)={l,0.5*l+0.5*a,0,h1};
Point(9)={0,0.5*l-0.5*a,0,h1};
Point(10)={b,0.5*l-0.5*a,0,h10};
Point(11)={b,0.5*l+0.5*a,0,h10};
Point(12)={0,0.5*l+0.5*a,0,h1};

Point(13)={0,0.5*l-25,0,h10};
Point(14)={l,0.5*l-25,0,h10};
Point(15)={l,0.5*l+25,0,h10};
Point(16)={0,0.5*l+25,0,h10};

Line(1)={1,2};
Line(2)={2,14};
Line(20)={14,5};
Line(3)={5,6};
Line(4)={6,7};
Line(5)={7,8};
Line(50)={8,15};
Line(6)={15,3};
Line(7)={3,4};
Line(70)={4,16};
Line(8)={16,12};
Line(9)={12,11};
Line(10)={11,10};
Line(11)={10,9};
Line(12)={9,13};
Line(13)={13,1};
Line(14)={14,13};
Line(15)={15,16};

Line Loop(1)={1,2,14,13};
Line Loop(2)={20,3,4,5,50,15,8,9,10,11,12,-14};
Line Loop(3)={6,7,70,-15};

Plane Surface(1)={1};
Plane Surface(2)={2};
Plane Surface(3)={3};

Physical Line(1)={1,2,20};  // fix both directions
Physical Line(2)={7};       // top edge
Physical Line(3)={70,8};    // left edge upper

Physical Surface(1)={1,2,3};


