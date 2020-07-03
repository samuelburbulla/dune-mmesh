lc = 0.1;
lcf = 0.01;

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

Physical Line(1) = {1:4};

// domain surface
Line Loop(1) = {1:4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};

// interface
Point(7) = {0.45, 0.5, 0, lcf};
Point(8) = {0.55, 0.5, 0, lcf};
Line(10) = {7, 8};
Line {10} In Surface{1};
Physical Line(10) = {10};

// interface boundary
Physical Point(1) = {7, 8};
