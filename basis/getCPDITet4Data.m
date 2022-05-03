function data = getCPDITet4Data(corners)
% Compute CPDI2 (linear tetrahedron) shape functions and first derivatives at
% particle "pid".
%
% Inputs:
%
% - corners: coordinates of the vertices of the tetrahedron
%
% VP Nguyen
% 20 April, 2016
% Monash University
% Clayton, Australia


% function and gradient weights

x1   = corners(1,1); y1   = corners(1,2); z1   = corners(1,3);
x2   = corners(2,1); y2   = corners(2,2); z2   = corners(2,3);
x3   = corners(3,1); y3   = corners(3,2); z3   = corners(3,3);
x4   = corners(4,1); y4   = corners(4,2); z4   = corners(4,3);

x12 = x1 - x2; x13 = x1 - x3; x14 = x1 - x4; x23 = x2 - x3; x24 = x2 - x4; x34 = x3 - x4;
x21 = -x12; x31=-x13;x41=-x14;x32=-x23;x42=-x24;x43=-x34;

y12 = y1 - y2; y13 = y1 - y3; y14 = y1 - y4; y23 = y2 - y3; y24 = y2 - y4; y34 = y3 - y4;
y21 = -y12; y31=-y13;y41=-y14;y32=-y23;y42=-y24;y43=-y34;

z12 = z1 - z2; z13 = z1 - z3; z14 = z1 - z4; z23 = z2 - z3; z24 = z2 - z4; z34 = z3 - z4;
z21 = -z12; z31=-z13;z41=-z14;z32=-z23;z42=-z24;z43=-z34;

det = 1/6;
V   = det*(x21*(y23*z34-y34*z23) + x32*(y34*z12-y12*z34) + x43*(y12*z23-y23*z12));
a1  = y42*z32 - y32*z42; b1 = x32*z42 - x42*z32; c1 = x42*y32 - x32*y42;
a2  = y31*z43 - y34*z13; b2 = x43*z31 - x13*z34; c2 = x31*y43 - x34*y13;
a3  = y24*z14 - y14*z24; b3 = x14*z24 - x24*z14; c3 = x24*y14 - x14*y24;
a4  = y13*z21 - y12*z31; b4 = x21*z13 - x31*z12; c4 = x13*y21 - x12*y31;

% V is a signed quantity!!!

%V = V;

wg    = zeros(4,3);
wf    = V*[1/4 1/4 1/4 1/4];


wg(1,:) = det*[a1 b1 c1];
wg(2,:) = det*[a2 b2 c2];
wg(3,:) = det*[a3 b3 c3];
wg(4,:) = det*[a4 b4 c4];

data.wf  = wf;
data.wg  = wg;
data.vol = V;


