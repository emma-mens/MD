function [Fij,Uij,rij]=force_calculation(N,r,Ls,rc2)

% Derivative of potential
   Fij = zeros(N,N,3);
 
   % Force b/n particle i and j
   % TODO (emazuh): optimize, account for r cutoff
   epsilon = 1.65*10^(-21); % LJ energy [J]
   sigma = 3.4*10^(-10); % LJ diameter [m]
   rij = zeros(N,N,3);
   rc = rc2^0.5;
   
   for i=1:N
       for j=1:N
           if i < j
               
               % set distance
               rij(i,j,:) = r(i)-r(j);
               rij(j,i,:) = - rij(i,j);
               
               % find magnitude and unit vector
               r_ij = abs(rij(i,j,:));
               
               % cut off radius
               if r_ij < rc
                   unit_vector = rij(i,j,:)/r_ij; 
                   s_over_r = sigma/r_ij;
                   
                   % calculate force (equal and opposite
                   Fij(i,j,:) = unit_vector*48*epsilon*((s_over_r)^12 - 0.5*(s_over_r)^6)/(r_ij^2);
                   Fij(j,i,:) = - Fij(i,j,:);
               end
           end
       end
   end
   
   Uij = zeros(N,N);
   
   

