function [Fij,Uij,rij]=force_calculation(N,r,Ls,rc2)
 
   % Force b/n particle i and j
   epsilon = 1.65*10^(-21); % LJ energy [J]
   sigma = 3.4*10^(-10); % LJ diameter [m]
   
   % calculate distance vectors and magnitudes b/n atom pairs
   rij  = permute(r - reshape(r', [1,3,N]), [1,3,2]);
   rij_mag = vecnorm(rij, 2, 3);


    % since I divide by the magnitude, I'll replace 0 magnitudes with 1
    % this is fine since the distance vector will still be zero for same atoms
    rij_mag(eye(N) == 1) = 1;

    s_r6 = (sigma./rij_mag).^6;
    s_r12 = s_r6.^2;

    Fij = (48 * epsilon * (s_r12 - 0.5*s_r6)./(rij_mag.^2)).*rij;


    Uij = 4*epsilon*(s_r12-s_r6);

    % because i hacked the distance matrix (0-> 1), i'll set self potentials to zero
    Uij(eye(N) == 1) = 0;


z


