function [r,v]=initialize(Ls,Ni,Ts)

N = Ni^3;


%*************************************************************************
% select initial velocities from Maxwell-Boltzmann distribution
%*************************************************************************

% From Wikipedia: https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
kB = 1.38065*10^(-23); % Boltzmann constant [J/K]
Ts = 1.171461; % desired (and initial) temperature in LJ units (Ts = kB*T/epsilon)
m = 39.948*1.660538921*10^(-27); % molecule mass [kg]
% a = (kB*T/m)^0.5;


% Using gaussian for the dist
sig = (kB*Ts/m)^0.5;
v = sig*randn(N, 3);

% Place atoms in a regular grid
r = ones(N,3);
particle_number = 1;
for i=1:Ls
    for j=1:Ls
        for k=1:Ls
            r(particle_number,:) = [i,j,k];
            particle_number = particle_number + 1;
        end
    end
end
            
