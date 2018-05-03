function [P,U]=P_and_U(N,Vs,rij,Fij,Uij,T_inst,P_LRC,E_LRC)    

% TODO (emazuh): Optimize this code
P_tensor = zeros(3,3);

for i=1:N
    for j=1:N
        if j>1
            P_tensor = P_tensor + rij(i,j)'*Fij(i,j);
        end
    end
end

P = (trace(P_tensor)/3 + P_LRC + N*kB*T_inst)/Vs;
U = E_LRC + mean2(Uij); %TODO: emazuh 