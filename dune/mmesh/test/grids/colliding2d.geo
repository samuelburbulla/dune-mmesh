lc = 0.1;
lcf = 0.03;

// domain corners
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

// domain outline
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// domain surface
Line Loop(1) = {1:4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};

// interface
Point(7) = {0.41, 0.25, 0, lcf};
Point(8) = {0.8, 0.25, 0, lcf};
Point(9) = {0.2, 0.75, 0, lcf};
Point(10) = {0.64, 0.75, 0, lcf};
Line(10) = {7, 8};
Line(11) = {9,10};
Line {10:11} In Surface{1};
Physical Line(10) = {10:11};
