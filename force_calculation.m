function [Fij,Uij,rij]=force_calculation(N,r,Ls,rc2)

   Uij = zeros(N,N);

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
   
%    rtest = reshape(r - reshape(r', [1,3,N]), [N,N,3]);
%    c = reshape(r', [1,3,N]);
%    r
%    r - c
%    c = reshape(r - c, N,N,3);
%    
%    c(2,1,:)
%    rij(2,1,:)
%    r(1,:) - c(1,:,2)
%    rtest(1,2,:)
%    t1 = rtest(1,1,:)
%    t2 = rij(1,1,:)

%    isequal(rij, reshape(rtest, [N,N,3]))
   
   

