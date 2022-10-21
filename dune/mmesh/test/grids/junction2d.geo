lc = 0.3;
lcf = 0.015;

// domain corners
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 0.5, 0, lcf};
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

// domain surface
Line Loop(1) = {1:6};
Plane Surface(1) = {1};
Physical Surface(1) = {1};

// interface
Point(7) = {0.6, 0.5, 0, lcf};
Point(8) = {0.75, 0.75, 0, lcf};
Line(10) = {3, 7};
Line(11) = {7, 6};
Line(12) = {7, 8};
Line {10:12} In Surface{1};
Physical Line(0) = {1:6};
Physical Line(1) = {10:12};
Physical Point(0) = {3, 6, 8};

// conflict point
Point(100) = {0.2, 0.52, 0, lcf};
Point(101) = {0.15, 0.51, 0, lcf};
Point{100:101} In Surface{1};
