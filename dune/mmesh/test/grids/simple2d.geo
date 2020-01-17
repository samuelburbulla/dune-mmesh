h = 2;

Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};
Point(5) = {1, 0.25, 0, h};
Point(6) = {0.5, 1, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 1};

Line(4) = {5, 3};
Line(5) = {3, 6};
Line(6) = {6, 5};

Line(7) = {6, 4};
Line(8) = {4, 1};
Line(9) = {1, 6};

Line(10) = {1, 5};
Line(11) = {5, 6};
Line(12) = {6, 1};

Line Loop(1) = {1:3};
Line Loop(2) = {4:6};
Line Loop(3) = {7:9};
Line Loop(4) = {10:12};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Physical Surface(1) = {1:4};
