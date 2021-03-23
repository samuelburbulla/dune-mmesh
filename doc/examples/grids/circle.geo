lc = 0.1;
lcf = 0.05;

// Domain
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 2, 0, lc};
Point(4) = {0, 2, 0, lc};

// Points of circle
Point(5)  = {0.5, 0.5, 0, lcf};
Point(6)  = {0.5, 0.7, 0, lcf};
Point(7)  = {0.3, 0.5, 0, lcf};
Point(8)  = {0.5, 0.3, 0, lcf};
Point(9)  = {0.7, 0.5, 0, lcf};

// Domain outline
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Interface outline
Ellipse(5) = {6, 5, 7};
Ellipse(6) = {7, 5, 8};
Ellipse(7) = {8, 5, 9};
Ellipse(8) = {9, 5, 6};

Curve Loop(1) = {5:8};
Curve Loop(2) = {1:4};

Plane Surface(0) = {2,1};
Plane Surface(1) = {1};
Physical Surface(0) = {0};
Physical Surface(1) = {1};

Physical Curve(10) = {5:8};
