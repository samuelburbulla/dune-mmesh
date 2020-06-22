lc = 0.05;
lcf = 0.05;

d = 1e-1;

// domain corners
Point(1) = { 0.0, 0.0, 0, lc};
Point(2) = { 1.0, 0.0, 0, lc};
Point(3) = { 1.0, 1.0, 0, lc};
Point(4) = { 0.5+0.5*d, 1.0, 0, lcf};
Point(5) = { 0.5-0.5*d, 1.0, 0, lcf};
Point(6) = { 0.0, 1.0, 0, lc};

// fracture corners
Point(10) = { 0.5+0.5*d, 0.5, 0, lcf};
Point(11) = { 0.5-0.5*d, 0.5, 0, lcf};

// domain outline
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

// fracture lines
Line(10) = {4, 10};
Line(11) = {10, 11};
Line(12) = {11, 5};

// domain surface
Line Loop(1) = {1:3, 10:12, 5:6};
Line Loop(2) = {-4, 10:12};

Plane Surface(0) = {1, 2};
Physical Surface(0) = {0};
Plane Surface(1) = {2};
Physical Surface(1) = {1};

// fracture domain boundary
Physical Line(10) = {10:12};
