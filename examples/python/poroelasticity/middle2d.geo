lc = 1;
lcf = 0.1;

// domain corners
Point(1) = {0, 0, 0, lc};
Point(2) = {80, 0, 0, lc};
Point(3) = {80, 80, 0, lc};
Point(4) = {0, 80, 0, lc};

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

// interface
Point(7) = {36, 40, 0, lcf};
Point(8) = {44, 40, 0, lcf};
Line(10) = {7, 8};
Line {10} In Surface{1};
Physical Line(10) = {10};

// interface boundary
Physical Point(1) = {7, 8};
