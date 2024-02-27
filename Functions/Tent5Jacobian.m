function J = Tent5Jacobian(Origin,Point1,Point2,Point3,Point4)

% The inputs to this function are the homogenous transformation matrices,
% the output is the Jacobian.

IDvec = [0;0;1];

R_00 = Origin(1:3,1:3);
R_01 = Point1(1:3,1:3);
R_02 = Point2(1:3,1:3);
R_03 = Point3(1:3,1:3);
%R_04 = Point4(1:3,1:3);

d_00 = Origin(1:3,4);
d_01 = Point1(1:3,4);
d_02 = Point2(1:3,4);
d_03 = Point3(1:3,4);
d_04 = Point4(1:3,4);

J = zeros(6,4);

% Column1
J(1:3,1) = cross((R_00*IDvec),(d_04-d_00));
J(4:6,1) = R_00*IDvec;

% Column2
J(1:3,2) = cross((R_01*IDvec),(d_04-d_01));
J(4:6,2) = R_01*IDvec;

% Column3
J(1:3,3) = cross((R_02*IDvec),(d_04-d_02));
J(4:6,3) = R_02*IDvec;

% Column3
J(1:3,4) = cross((R_03*IDvec),(d_04-d_03));
J(4:6,4) = R_03*IDvec;



end