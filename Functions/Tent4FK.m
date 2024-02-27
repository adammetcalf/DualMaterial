function [OriginDH_x, OriginRot_z, OriginDH_y, OriginRot_Negz,...
    Point1DH, Point2DH, Point3DH, Point4DH] = Tent4FK(alpha0,beta0)


OriginDH_x = [1 0 0 0
              0 cos(alpha0) -sin(alpha0) 0
              0 sin(alpha0) cos(alpha0) 0
              0 0 0 1];

OriginRot_z = [cos(-pi/2) -sin(-pi/2) 0 0
               sin(-pi/2) cos(-pi/2) 0 0
               0 0 1 0
               0 0 0 1];

OriginDH_y = [1 0 0 0
              0 cos(beta0) -sin(beta0) 0
              0 sin(beta0) cos(beta0) 0
              0 0 0 1];

OriginRot_Negz = [cos(pi/2) -sin(pi/2) 0 0
               sin(pi/2) cos(pi/2) 0 0
               0 0 1 0
               0 0 0 1];

theta1 = deg2rad(0);
alpha1 = deg2rad(0);
Beta1 = deg2rad(0);

Point1DH = [cos(0) -sin(0)*cos(alpha1) sin(0)*sin(alpha1) 0
            sin(0) cos(0)*cos(alpha1) -cos(0)*sin(alpha1) 0
            0 sin(alpha1) cos(alpha1) 0.01
            0 0 0 1];

theta2 = deg2rad(0);
alpha2 = deg2rad(0);

Point2DH = [cos(theta2) -sin(theta2)*cos(alpha2) sin(theta2)*sin(alpha2) 0
            sin(theta2) cos(theta2)*cos(alpha2) -cos(theta2)*sin(alpha2) 0
            0 sin(alpha2) cos(alpha2) 0.01
            0 0 0 1];

theta3 = deg2rad(0);
alpha3 = deg2rad(0);

Point3DH = [cos(theta3) -sin(theta3)*cos(alpha3) sin(theta3)*sin(alpha3) 0
            sin(theta3) cos(theta3)*cos(alpha3) -cos(theta3)*sin(alpha3) 0
            0 sin(alpha3) cos(alpha3) 0.01
            0 0 0 1];

theta4 = deg2rad(0);
alpha4 = deg2rad(0);

Point4DH = [cos(theta4) -sin(theta4)*cos(alpha4) sin(theta4)*sin(alpha4) 0
            sin(theta4) cos(theta4)*cos(alpha4) -cos(theta4)*sin(alpha4) 0
            0 sin(alpha4) cos(alpha4) 0.01
            0 0 0 1];

end

