lc = 0.1;
lcf = 0.1;

// domain corners
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
Point(5) = {0, 0, 1, lc};
Point(6) = {1, 0, 1, lc};
Point(7) = {1, 1, 1, lc};
Point(8) = {0, 1, 1, lc};
Point(9) = {0.3, 0.55, 0.2, lc};
Point(10) = {0, 0.5, 0, lcf};
Point(11) = {1, 0.5, 0, lcf};
Point(12) = {1, 0.5, 1, lcf};
Point(13) = {0, 0.5, 1, lcf};

// domain outline
Line(1) = {8, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 7};
//+
Line(4) = {7, 8};
//+
Line(5) = {5, 1};
//+
Line(6) = {1, 2};
//+
Line(7) = {2, 6};
//+
Line(8) = {6, 5};
//+
Line(9) = {5, 13};
//+
Line(10) = {13, 8};
//+
Line(11) = {1, 10};
//+
Line(12) = {10, 4};
//+
Line(13) = {2, 11};
//+
Line(14) = {11, 3};
//+
Line(15) = {6, 12};
//+
Line(16) = {12, 7};
//+
Line(17) = {12, 13};
//+
Line(18) = {13, 10};
//+
Line(19) = {10, 11};
//+
Line(20) = {11, 12};
//+
Curve Loop(1) = {18, 19, 20, 17};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 2, 3, 4};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, 5, 6, 7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {9, 18, -11, -5};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {7, 15, -20, -13};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {8, 9, -17, -15};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {6, 13, -19, -11};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {14, -2, -12, 19};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {20, 16, -3, -14};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {17, 10, -4, -16};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {18, 12, -1, -10};
//+
Plane Surface(11) = {11};
//+
Surface Loop(1) = {11, 8, 9, 10, 2, 1};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {4, 6, 3, 7, 5, 1};
//+
Volume(2) = {2};

Point{9} In Volume{1};

Physical Volume(1) = {1,2};
Physical Surface(1) = {1};
