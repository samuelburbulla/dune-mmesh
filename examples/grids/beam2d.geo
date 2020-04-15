lc = 0.05;

// domain corners
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 0.15, 0, lc};
Point(4) = {0, 0.15, 0, lc};

// domain outline
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Physical Line(1) = {1:4};

// domain surface
Line Loop(1) = {1:4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};
