function [Fij,Uij,rij]=force_calculation(N,r,Ls,rc2)
   
   % calculate distance vectors and magnitudes b/n atom pairs
   rij_init  = permute(r - reshape(r', [1,3,N]), [1,3,2]);
   rp = rij_init + Ls;
   rm = rij_init - Ls;
   
   % update rij for closes mirror of other atoms in the box
   rij = rij_init.*(rij_init <= 0.5*Ls) + rm.*(rij_init > 0.5*Ls);
   rij = rij.*(rij >= -0.5*Ls) + rp.*(rij < -0.5*Ls);
   
   rij_mag = vecnorm(rij, 2, 3);

   % since I divide by the magnitude, I'll replace 0 elf distance magnitudes with 1
   % this is fine since s_r12 - s_r6 will still be zero for same atoms
   rij_mag(eye(N) == 1) = 1;

   s_r6 = (1./rij_mag).^6;
   s_r12 = s_r6.^2;

   Fij = (48*(s_r12 - 0.5*s_r6)./(rij_mag.^2)).*rij;

   Uij = 4*(s_r12-s_r6);
