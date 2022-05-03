a = 4;
b = 1;
t = 1;

h  = 0.09;
h1 = 0.1*h;
h2 = 0.4*h;

// -------------------
// corner points 
// -------------------

Point(1) = {0,0,0,h};
Point(2) = {a,0,0,h};
Point(3) = {a,b,0,h};
Point(4) = {0.,b,0,h};


// -------------------
// outer lines
// -------------------

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


Line Loop(100) = {1,2,3,4};


// -------------------
//  Surfaces
// -------------------

Plane Surface(100) = {100}; 


//Recombine Surface{100};
//Recombine Surface{200};

// lower, right, upper, left edges
// opposite edges = same direction

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};

Extrude {0, 0, t} {
  Surface{100};
}



// ----------------------
// Physical quantities
// ----------------------

Physical Volume(1) = {1};


