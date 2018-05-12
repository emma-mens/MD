function [P,U]=P_and_U(N,Vs,rij,Fij,Uij,T_inst,P_LRC,E_LRC)    

kB = 1.38065*10^(-23); % Boltzmann constant [J/K]

% TODO (emazuh): Optimize this code
Pij = zeros(3,3);

% optimize later

for i=1:3
    for j=1:3
        for k=1:N
            for l=k:N
                Pij(i,j) = Pij(i,j) + rij(k,l,i)*Fij(k,l,j);
            end
        end
    end
end

% for k=1:N
%     for l=1:N
%         if l > k
%             % P_tensor = P_tensor + rij(i,j)'*Fij(i,j);
%             Pij(i,j) = 
%         end
%     end
% end

P = trace(Pij/Vs)/3 + P_LRC + N*kB*T_inst/Vs
U = E_LRC + mean2(Uij) %TODO: emazuh 
pause on