%*************************************************************************
%*************************************************************************
%*************************************************************************
%
%                   2.37 Fundamentals of Nanoengineering
%                   Molecular Dynamics Project
%                   Originally written by Michael S. H. Boutilier
%                   Modified by N. G. Hadjiconstantinou 
%                   Aug. 1, 2014
%                   Submitted by: Emmanuel Azuh
%
%*************************************************************************
%*************************************************************************
%*************************************************************************

clear all
clc;

%*************************************************************************
%*************************************************************************
%*************************************************************************
%
% setup simulation parameters
%
%*************************************************************************
%*************************************************************************
%*************************************************************************


%*************************************************************************
% input parameters
%*************************************************************************
Ni = 8; %8 number of atoms per side in original cubic configuration 
Ts = 2; % 1.171461; % desired (and initial) temperature in LJ units (Ts = kB*T/epsilon)
ns = 0.9; % 0.903992; % LJ number density ns = n*sigma^3
epsilon = 1.65*10^(-21); % LJ energy [J]
sigma = 3.4*10^(-10); % LJ diameter [m]
m = 39.948*1.660538921*10^(-27); % molecule mass [kg]

kB = 1.38065*10^(-23); % Boltzmann constant [J/K]

%*************************************************************************
% solver parameters
%*************************************************************************

rc = 3; % cut-off radius, in LJ units, i.e., number of sigmas
STEPS = 1*10^4; %4 total number of time steps
STEPS_thermostat = 3*10^3; % number of steps to leave thermostat on until
STEPS_equilib = 4*10^3; % number of steps before starting to average
alpha = 0.01; % constant temperature constraint method relaxation parameter (optional)
dt = 0.005; % time step, in LJ units

%*************************************************************************
% calculated system parameters
%*************************************************************************

N = Ni^3; % number of atoms 
Vs = N/ns; % volume of domain, in LJ units, Vs = V/sigma^3
Ls = Vs^(1/3); % length of the domain in each direction, in LJ units
T_set = Ts*epsilon/kB; % set temperature, in Kelvin
P_LRC = 32/9*pi*ns^2*rc^(-9) - 16/3*pi*ns^2*rc^(-3); % long-range P correction
E_LRC = 8/9*pi*ns*rc^(-9)-8/3*pi*ns*rc^(-3); % long-range U correction (per particle)
rc2 = rc^2; % square of cut-off radius


%*************************************************************************
% initialize variables
%*************************************************************************

% r = zeros(N,3); % atom position vectors at current time---3 components
rij = zeros(N,N,3); % separation (matrix) for atom i and j at current time---3 components
% v = zeros(N,3); % atom velocity vectors---3 components 
F = zeros(N,3); % total force on each atom at current time---3 components
Res = zeros(STEPS,4); % Result Array (for all timesteps)  

%*************************************************************************
%*************************************************************************
%*************************************************************************
%
% define initial configuration (r,v for all N atoms)
% Recommend to use input variable Ni to set up a cubic lattice for positions
% and a Maxwell-Boltzmann distribution at temperature Ts for the velocities
%
%*************************************************************************
%*************************************************************************
%*************************************************************************
[r,v]=initialize(Ls,Ni,Ts);
init_pos = r;

% %*************************************************************************
% % Remove any center of mass motion
% %************************************************************************
c = sum(v)/N;
v = v-repmat(c,N,1);

%*************************************************************************
%*************************************************************************
%*************************************************************************
%
% Solve equations of motion as a loop over STEPS timesteps
%
%*************************************************************************
%*************************************************************************
%*************************************************************************

for t = 1:1:STEPS   
    
%*************************************************************************
% calculate forces 
%*************************************************************************

%*************************************************************************
% routine "force_calculation" returns matrix Fij of size NxNx3 
% containing the force (3 components) between particle i and j and matrix
% Uij of size NxN containing the corresponding energy of interaction
% *************************************************************************
[Fij,Uij,rij]= force_calculation(N,r,Ls,rc2);

    F = squeeze(sum(Fij, 2)); % sum along second dimension

%*************************************************************************
% calculate new atom positions, velocities based on forces just calculates
% scale velocities such that T=Ts if t < STEPS_thermostat
% apply periodic boundary conditions for atoms that "step out" of the box
%*************************************************************************
[rnew,vnew,T_inst]= take_one_step(N,r,v,F,t,dt,Ls,STEPS_thermostat,alpha,Ts);
    
    r = rnew;
    v = vnew;
    
%*************************************************************************
% compute current pressure and potential energy
%*************************************************************************
    [P,U]= P_and_U(N,Vs,rij,Fij,Uij,T_inst,P_LRC,E_LRC);
    
    Res(t,:) = [t T_inst P U];

    if mod(t,100)==0
        [t T_inst P U] % write to screen        
    end
    
end

%*************************************************************************
%*************************************************************************
%*************************************************************************
%
% post-processing
%
%*************************************************************************
%*************************************************************************
%*************************************************************************



figure(1), hold on
subplot(3,1,1), hold on
plot([1:size(Res,1)],Res(:,2),'k'), hold on
subplot(3,1,2), hold on
plot([1:size(Res,1)],Res(:,3),'k'), hold on
subplot(3,1,3), hold on
plot([1:size(Res,1)],Res(:,4),'k'), hold on

N_start = STEPS_equilib;
[mean(Res(N_start:size(Res,1),2)) mean(Res(N_start:size(Res,1),3)) mean(Res(N_start:size(Res,1),4))]
 
% D = mean(norm(x)^2)/(6*Ts)

% Results:
% rho .5, T 5 -> [T, P, U ] 4.9598    4.6547   -2.3422
% rho .9  T 2 -> [T, P, U]  1.9124    8.7789   -5.0390
% rho .8  T 4 -> [T, P, U]  3.9715   12.1735   -3.4477

