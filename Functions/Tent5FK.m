function[Origin, Point1, Point2, Point3, Point4,J] = Tent5FK(theta0,alpha0,alpha1,alpha2,alpha3,alpha4,StartTransform,l1,l2,l3,l4)

OriginDH = [cos(theta0) -sin(theta0)*cos(alpha0) sin(theta0)*sin(alpha0) 0
            sin(theta0) cos(theta0)*cos(alpha0) -cos(theta0)*sin(alpha0) 0
            0 sin(alpha0) cos(alpha0) 0
            0 0 0 1];

Point1DH = [1 0 0 0
            0 cos(alpha1) -sin(alpha1) 0
            0 sin(alpha1) cos(alpha1) l1
            0 0 0 1];

Point2DH = [1 0 0 0
            0 cos(alpha2) -sin(alpha2) 0
            0 sin(alpha2) cos(alpha2) l2
            0 0 0 1];

Point3DH = [1 0 0 0
            0 cos(alpha3) -sin(alpha3) 0
            0 sin(alpha3) cos(alpha3) l3
            0 0 0 1];

Point4DH = [1 0 0 0
            0 cos(alpha4) -sin(alpha4) 0
            0 sin(alpha4) cos(alpha4) l4
            0 0 0 1];


% Evaluate Forward Kinematics (and homogeneous transformation matrices)
Origin = StartTransform*OriginDH;
Point1 = Origin*Point1DH;
Point2 = Point1*Point2DH;
Point3 = Point2*Point3DH;
Point4 = Point3*Point4DH;

% Evaluate Jacobian
J = Tent5Jacobian(Origin,Point1,Point2,Point3,Point4);
end