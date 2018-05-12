function [Fij,Uij,rij]=force_calculation(N,r,Ls,rc2)

   Uij = zeros(N,N);

   % Derivative of potential
   Fij = zeros(N,N,3);
 
   % Force b/n particle i and j
   epsilon = 1.65*10^(-21); % LJ energy [J]
   sigma = 3.4*10^(-10); % LJ diameter [m]
   rc = rc2^0.5;
   
   % calculate distance b/n atom pairs
   rij  = permute(r - reshape(r', [1,3,N]), [1,3,2]);
   
   % TODO (emazuh): optimize force calculation, account for r cutoff
   for i=1:N
       for j=1:N
           if i < j
               
               % set distance
               rij(i,j,:) = r(i,:)-r(j,:);
               rij(j,i,:) = - rij(i,j,:);

               % find magnitude
               r_ij = norm(r(i,:)-r(j,:));
               
               % cut off radius
               if r_ij < rc
                   s_over_r = sigma/r_ij;
                   
                   % calculate force (equal and opposite)
                   Fij(i,j,:) = rij(i,j,:)*48*epsilon*((s_over_r)^12 - 0.5*(s_over_r)^6)/(r_ij^2);
                   Fij(j,i,:) = - Fij(i,j,:);
                   
                   % LJ potential
                   Uij(i,j) = 4*epsilon*(s_over_r^12 - s_over_r^6);
                   Uij(j, i) = Uij(i, j);
               end
           end
       end
   end
  


