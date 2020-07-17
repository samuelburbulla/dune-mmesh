lc = 0.1;
lcf = lc;

// Domain
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 2, 0, lc};
Point(4) = {0, 2, 0, lc};

// Points of circle
Point(6)  = {0.3, 0.3, 0, lcf};
Point(7)  = {0.7, 0.3, 0, lcf};
Point(8)  = {0.7, 0.7, 0, lcf};
Point(9)  = {0.3, 0.7, 0, lcf};

// Domain outline
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Interface outline
Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 6};

Curve Loop(1) = {5:8};
Curve Loop(2) = {1:4};

Plane Surface(1) = {1};
Plane Surface(2) = {2,1};
Physical Surface(1) = {1};
Physical Surface(2) = {2};

Physical Curve(10) = {5:8};
