lc = 0.1;
lcf = 0.1;

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
Point(7) = {0.1, 0.5, 0, lcf};
Point(8) = {0.5, 0.3, 0, lcf};
Point(9) = {0.3, 0.65, 0, lcf};
Line(11) = {6,7};
Line(12) = {8,9};
Line {11:12} In Surface{1};
Physical Line(1) = {11:12};
