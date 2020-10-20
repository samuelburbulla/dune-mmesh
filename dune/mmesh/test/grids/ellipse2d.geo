//characteristic length
lc = 5e-2;

//domain corners, domain = [0,1] x [0,1]
Point(1) = {-0.5, -0.5, 0, lc};
Point(2) = {0.5, -0.5, 0, lc};
Point(3) = {0.5, 0.5, 0, lc};
Point(4) = {-0.5, 0.5, 0, lc};

Point(5) = {0, 0, 0, lc}; //ellipse center
Point(6) = {0.3, 0, 0, lc};
Point(7) = {-0.3, 0, 0, lc};
Point(8) = {0, 0.15, 0, lc};
Point(9) = {0, -0.15, 0, lc};


//domain outline
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//interface
Ellipse(5) = {6, 5, 6, 8}; //ellipse arc with start point 6, center 5, point on major axis 6, end point 8
Ellipse(6) = {8, 5, 6, 7};
Ellipse(7) = {7, 5, 6, 9};
Ellipse(8) = {9, 5, 6, 6};

//curve loops
Curve Loop(1) = {1:4}; //domain boundary

//surfaces
Plane Surface(1) = {1};
Physical Surface(1) = {1};

Curve{5:8} In Surface{1};
Physical Curve(2) = {5:8};
