lc = 0.1;
lcf = 0.1;
lcg = 0.1;

// domain corners
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 0.5, 0, lcg};
Point(4) = {1, 1, 0, lc};
Point(5) = {0, 1, 0, lc};
Point(6) = {0, 0.5, 0, lcf};

// domain outline
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

// domain boundary
Physical Line(1) = {1:6};

// domain surface
Line Loop(1) = {1:6};
Plane Surface(1) = {1};
Physical Surface(1) = {1};

// interface
Line(12) = {3, 6};
Line {12} In Surface{1};
Physical Line(2) = {12};

// interface boundary
Physical Point(1) = {3, 6};
