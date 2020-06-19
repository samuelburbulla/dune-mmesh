lc = 0.1;
lcf = 0.01;

// Domain
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

// Points of circle
Point(5)  = {0.5, 0.5, 0, lcf};
Point(6)  = {0.5, 0.7, 0, lcf};
Point(7)  = {0.3, 0.5, 0, lcf};
Point(8)  = {0.5, 0.3, 0, lcf};
Point(9)  = {0.7, 0.3, 0, lcf};
Point(10) = {0.7, 0.9, 0, lcf};
Point(11) = {0.5, 0.9, 0, lcf};

// Domain outline
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Interface outline
Ellipse(5) = {6, 5, 6, 7};
Ellipse(6) = {7, 5, 7, 8};
Line(7) = {8, 9};
Ellipse(8) = {6, 11, 6, 10};

Curve Loop(1) = {1:4};

Plane Surface(1) = {1};
Physical Surface(1) = {1};

Physical Curve(10) = {5:8};
Curve{5:8} In Surface{1};
