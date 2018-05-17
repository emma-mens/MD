function [P,U]=P_and_U(N,Vs,rij,Fij,Uij,T_inst,P_LRC,E_LRC)    

% since we only care about the xx,yy and zz components, I only calculate those
% Divide by 2 to account for double counting
Pij = rij.*Fij;
P_axes = sum(Pij(:))/3/Vs/2;
P = P_axes + N*T_inst/Vs; 

% Energy is double counted since i calculate for all pairs so I divide by 2
U = sum(Uij(:))/N/2; % LJ
