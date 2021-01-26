//characteristic length
lc = 9e-2;
lc2 = 9e-2;

// ------------------- domain -------------------
//domain = [-0.5,0.5] x [-0.5,0.5] x [-0.5,0.5]

//domain corners
Point(101) = {-0.5, -0.5, -0.5, lc2};
Point(102) = { 0.5, -0.5, -0.5, lc2};
Point(103) = { 0.5,  0.5, -0.5, lc2};
Point(104) = {-0.5,  0.5, -0.5, lc2};
Point(105) = {-0.5, -0.5,  0.5, lc2};
Point(106) = { 0.5, -0.5,  0.5, lc2};
Point(107) = { 0.5,  0.5,  0.5, lc2};
Point(108) = {-0.5,  0.5,  0.5, lc2};

//domain edges
Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 104};
Line(104) = {104, 101};
Line(105) = {105, 106};
Line(106) = {106, 107};
Line(107) = {107, 108};
Line(108) = {108, 105};
Line(109) = {101, 105};
Line(110) = {102, 106};
Line(111) = {103, 107};
Line(112) = {104, 108};

//domain surfaces
Line Loop(113) = {101:104};
Plane Surface(114) = {113};
Line Loop(115) = {105:108};
Plane Surface(116) = {115};
Line Loop(117) = {102, 111, -106, -110};
Plane Surface(118) = {117};
Line Loop(119) = {103, 112, -107, -111};
Plane Surface(120) = {119};
Line Loop(121) = {104, 109, -108, -112};
Plane Surface(122) = {121};
Line Loop(123) = {101, 110, -105, -109};
Plane Surface(124) = {123};

Surface Loop(125) = {114, 116, 118, 120, 122, 124};
Volume(126) = {125};
Physical Volume(1) = {126};


// ------------------- sphere -------------------

//major axes
a = 0.35;
b = 0.15;
c = 0.25;

//points for interface (ellisoid)
Point(1) = {0, 0, 0, lc}; //center
Point(2) = {a, 0, 0, lc};
Point(3) = {0, b, 0, lc};
Point(4) = {-a, 0, 0, lc};
Point(5) = {0, -b, 0, lc};
Point(6) = {0, 0, c, lc};
Point(7) = {0, 0, -c, lc};

//interface outline
//ellipsis in xy-plane
Ellipse(1) = {2, 1, 2, 3}; //ellipse arc with start point 2, center 1, point on major axis 2, end point 3
Ellipse(2) = {3, 1, 2, 4};
Ellipse(3) = {4, 1, 2, 5};
Ellipse(4) = {5, 1, 2, 2};

//ellipsis in xz-plane
Ellipse(5) = {6, 1, 2, 2};
Ellipse(6) = {2, 1, 2, 7};
Ellipse(7) = {7, 1, 2, 4};
Ellipse(8) = {4, 1, 2, 6};

//ellipsis in yz-plane
Ellipse(9) =  {6, 1, 6, 5};
Ellipse(10) = {5, 1, 6, 7};
Ellipse(11) = {7, 1, 6, 3};
Ellipse(12) = {3, 1, 6, 6};

//sphere surfaces
Line Loop(13) = {1, 12, 5};
Surface(14) = {13};
Line Loop(15) = {12, -8, -2};
Surface(16) = {15};
Line Loop(17) = {8, 9, -3};
Surface(18) = {17};
Line Loop(19) = {9, 4, -5};
Surface(20) = {19};
Line Loop(21) = {1, -11, -6};
Surface(22) = {21};
Line Loop(23) = {10, -6, -4};
Surface(24) = {23};
Line Loop(25) = {10, 7, 3};
Surface(26) = {25};
Line Loop(27) = {11, 2, -7};
Surface(28) = {27};

Surface Loop(29) = {16, 14, 22, 28, 26, 24, 20, 18};
Surface{16, 14, 22, 28, 26, 24, 20, 18} In Volume{126};
Physical Surface(10) = {16, 14, 22, 28, 26, 24, 20, 18};
