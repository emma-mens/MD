function [Fij,Uij,rij]=force_calculation(N,r,Ls,rc2)
 
   % epsilon = 1.65*10^(-21); % LJ energy [J]
   %sigma = 3.4*10^(-10); % LJ diameter [m]
   
   % calculate distance vectors and magnitudes b/n atom pairs
   rij_init  = permute(r - reshape(r', [1,3,N]), [1,3,2]);
   rp = permute(r - reshape((r+Ls)', [1,3,N]), [1,3,2]);
   rm = permute(r - reshape((r-Ls)', [1,3,N]), [1,3,2]);
   rfinal = rij_init.*(abs(rij_init) < abs(rp)) + rp.*(abs(rp) < abs(rij_init));
   rij = rfinal.*(abs(rfinal) < abs(rm)) + rm.*(abs(rm) < abs(rfinal));
   
   % set self distances to zero and distances between corner atoms to
   % rij (without any offset since an offset will cause overlapping)
   rij(rij == 0) = rij_init(rij == 0);
   rij(eye(N) == 1) = 0;
   
   rij_mag = vecnorm(rij, 2, 3)

   % rx_mag = 

   % since I divide by the magnitude, I'll replace 0 magnitudes with 1
   % this is fine since s_r12 - s_r6 will still be zero for same atoms
   rij_mag(eye(N) == 1) = 1;

   s_r6 = (1./rij_mag).^6;
   s_r12 = s_r6.^2;

   Fij = (48*(s_r12 - 0.5*s_r6)./(rij_mag.^2)).*rij;

   Uij = 4*(s_r12-s_r6);



%    rx = r(:,1); ry = r(:,2); rz = r(:,3);
%    rij = zeros(N,N,3);
%    
%    % rij between mirros
%    % For x: do x + Ls, x - Ls,
%    rxplus = rx + Ls; rxminus = rx - Ls; ryplus = ry + Ls; ryminus = ry - Ls;
%    rzplus = rz + Ls; rzminus = rz - Ls;
% 
%    rijx = permute(rx - reshape(rx', [1,1,N]), [1,3,2]);
%    rijxp = permute(rx - reshape(rxplus', [1,1,N]), [1,3,2]);
%    rijxm = permute(rx - reshape(rxminus', [1,1,N]), [1,3,2]);
%    
%    rijy = permute(ry - reshape(ry', [1,1,N]), [1,3,2]);
%    rijyp = permute(ry - reshape(ryplus', [1,1,N]), [1,3,2]);
%    rijym = permute(ry - reshape(ryminus', [1,1,N]), [1,3,2]);
%    
%    % select images closest to any given atom
%    tmpx = rijx.*(abs(rijx) < abs(rijxp)) + rijxp.*(abs(rijxp) < abs(rijx));
%    tmpx = tmpx.*(abs(tmpx) < abs(rijxm)) + rijxm.*(abs(rijxm) < abs(tmpx))
%    
%    tmpy = rijy.*(abs(rijy) < abs(rijyp)) + rijyp.*(abs(rijyp) < abs(rijy));
%    tmpy = tmpy.*(abs(tmpy) < abs(rijym)) + rijym.*(abs(rijym) < abs(tmpy))
