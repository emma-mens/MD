function [rnew,vnew,T_inst]=take_one_step(N,r,v,F,t,dt,Ls,STEPS_thermostat,alpha,Ts)

m = 39.948*1.660538921*10^(-27); % molecule mass [kg]
kB = 1.38065*10^(-23); % Boltzmann constant [J/K]
sigma = 3.4*10^(-10); % LJ diameter [m]
epsilon = 1.65*10^(-21); % LJ energy [J]
rc = 3;
rc2 = rc^2;
% Ls = Ls*sigma; % TODO (emazuh): change this if you use LJ Ls in force_calculation

% calculate acceleration
a = F/m;


% TODO (emazuh): velocity verley (maybe switch with leapfrog later?)

% update positions
rnew = r + v*dt + 0.5*a*dt^2;

% apply boundary condition on rnew
rnew(rnew < 0 | rnew > Ls) = mod(rnew(rnew < 0 | rnew > Ls), Ls);

[Fij, ~, ~] = force_calculation(N,rnew,Ls, rc2);
F_next = squeeze(sum(Fij, 2));

% update velocities
vnew = v + 0.5*(a + F_next/m)*dt;

% TODO (emazuh): calculate temperature (using thermostate with large Q)
%V = (Ls*sigma)^3;
%rho = m*N/V;
% u = m*sum(vnew, 1)/rho/V; % flow velocity
u = mean(vnew, 1);

% Ts = 1.171461; % desired (and initial) temperature in LJ units (Ts = kB*T/epsilon)
% T_set = Ts*epsilon/kB; % set temperature, in Kelvin
% TODO (emazuh): apply thermostat (also check if t < STEPS_thermostat
% d(zeta)/dt = 6NkBT/Q *(T_set - T_inst)
% ma = f - zeta*p

% p = m*norm(vnew).^2-norm(u)^2;

sig = (kB*Ts/m)^0.5;


% Use crude method for now. TODO (emazuh): implement Nose-Hoover
if (mod(t, 200)) == 0 && (t < STEPS_thermostat)
    vnew = sig*randn(N, 3);
end

T_inst = m*sum(vecnorm(vnew, 2, 2).^2-norm(u)^2, 1)/(3*N*kB);
T_inst = T_inst*kB/epsilon; % LJ    
