function[Mass,I,k,Length] = MaterialProperties(radius,rho,DH,E)

% Link Lengths (m)
Length = DH(3,4);


% Link Volumes (m^3)
Vol = Length*pi*(radius^2);


% Link Masses (consider as a point mass, with the mass of each link acting
% at the endpoint - TODO: Adjust this to act at the COM of each link.)
Mass = rho*Vol;


% Define Inertia (kg m^-2) 
%%% TO DO: Inertia evaluated about correct position?
%%% SOLID CYCLNDER ABOUT ENDPOINT DIAMETER
%%% I = 1/4*m*R^2 + 1/3*m*L^2

I = (1/4)*Mass*radius^2 + (1/3)*Mass*Length;


% Define Joint Stiffness (Nm rad^-1)
% k_m = EI L^-1

% NOTE, stiffness evaluated for Link01 is applied at joint0 etc

k = (E*I)/Length;


end