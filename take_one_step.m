function [rnew,vnew,T_inst]=take_one_step(N,r,v,F,t,dt,Ls,STEPS_thermostat,alpha,Ts)

    m = 1; 

    rc = 3;
    rc2 = rc^2;

    % calculate acceleration
    a = F/m;

    % update positions
    v = v + 0.5*a*dt;
%     rnew = r + vhalf*dt;
    rnew = r + v*dt + 0.5*a*dt^2;

    % apply boundary condition on rnew
    rnew(rnew < 0 | rnew > Ls) = mod(rnew(rnew < 0 | rnew > Ls), Ls);
    for i=1:N
        for j=1:3
            if rnew(i,j) < 0 | rnew(i,j) > Ls
                i
                j
                rnew(i,j)
                error('moji')
            end
        end
    end
%     [Fij, ~, ~] = force_calculation(N,rnew,Ls, rc2);
%     F_next = squeeze(sum(Fij, 2));

    % update velocities
%     anew = F_next/m;
%     vnew = v + 0.5*(a + anew)*dt;
    
    

    % assert(rnew(rnew >= 0) & rnew(rnew <= Ls), 'atoms are outside the box');

    u = mean(v, 1);

    T_inst = m*sum(vecnorm(v, 2, 2).^2-norm(u)^2, 1)/(3*N);

    if t < STEPS_thermostat
       v = v*sqrt(Ts/T_inst); 
    end
    vnew = v + 0.5*a*dt;
