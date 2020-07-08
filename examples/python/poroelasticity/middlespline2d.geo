lc = 5;
lcf = 0.1;

// domain corners
Point(1) = {0, 0, 0, lc};
Point(2) = {100, 0, 0, lc};
Point(3) = {100, 100, 0, lc};
Point(4) = {0, 100, 0, lc};

// domain outline
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Physical Line(1) = {1:4};

// fracture
Point(7) = {45, 50, 0, lcf};
Point(8) = {55, 50, 0, lcf};
Point(9) = {50, 49.5, 0, lcf};
Point(10) = {50, 50.5, 0, lcf};
BSpline(10) = {7, 9, 8};
BSpline(11) = {8, 10, 7};


// domain surface
Line Loop(1) = {1:4};
Line Loop(2) = {10:11};
Plane Surface(1) = {1,2};
Plane Surface(2) = {2};
Physical Surface(0) = {1};
Physical Surface(1) = {2};
