function [Fij,Uij,rij]=force_calculation(N,r,Ls,rc2)

% Derivative of potential
   size(r)
   Fij = zeros(N,N,3);
   
   %% 
   % Force b/n particle i and j
   fij = 48*(rij^(-14) - 0.5*);
   
   %%
   Uij = zeros(N,N);
   rij = zeros(N,N,3);
   

