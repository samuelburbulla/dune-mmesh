lc = 0.5;
lcf = 0.1;
d = 1.0;

// domain corners
Point(1) = {0,   0, 0, lc};
Point(2) = {100,  0, 0, lc};
Point(3) = {100, 100, 0, lc};
Point(4) = {0,  100, 0, lc};

// domain outline
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Physical Line(1) = {1:4};

// fracture
Point(7)  = {45, 50-0.5*d, 0, lcf};
Point(8)  = {55, 50-0.5*d, 0, lcf};
Point(9)  = {45, 50+0.5*d, 0, lcf};
Point(10) = {55, 50+0.5*d, 0, lcf};
Line(10) = {7, 8};
Line(11) = {8, 10};
Line(12) = {10, 9};
Line(13) = {9, 7};


// domain surface
Line Loop(1) = {1:4};
Line Loop(2) = {10:13};
Plane Surface(1) = {1,2};
Plane Surface(2) = {2};
Physical Surface(0) = {1};
Physical Surface(1) = {2};
