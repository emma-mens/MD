function [P,U]=P_and_U(N,Vs,rij,Fij,Uij,T_inst,P_LRC,E_LRC)    

kB = 1.38065*10^(-23); % Boltzmann constant [J/K]
epsilon = 1.65*10^(-21); % LJ energy [J]
sigma = 3.4*10^(-10); % LJ diameter [m]

% since we only care about the xx,yy and zz components, I only calculate those
% Divide by 2 to account for double counting
Pij = rij.*Fij;
P_axes = sum(Pij(:))/3/Vs/2;
P = P_axes + N*kB*T_inst/Vs;
P = P*sigma^3/epsilon % LJ

U = sum(Uij(:))/epsilon % LJ
z
pause on