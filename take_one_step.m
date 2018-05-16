function [rnew,vnew,T_inst]=take_one_step(N,r,v,F,t,dt,Ls,STEPS_thermostat,alpha,Ts)

    m = 1; 

    rc = 3;
    rc2 = rc^2;

    % calculate acceleration
    a = F/m;

    % update positions
    vhalf = v + 0.5*a*dt;
    rnew = r + vhalf*dt;

    % apply boundary condition on rnew
    rnew(rnew < 0 | rnew > Ls) = mod(rnew(rnew < 0 | rnew > Ls), Ls);

    [Fij, ~, ~] = force_calculation(N,rnew,Ls, rc2);
    F_next = squeeze(sum(Fij, 2));
    z
    % update velocities
    anew = F_next/m;
    vnew = v + 0.5*(a + anew)*dt;

    % assert(rnew(rnew >= 0) & rnew(rnew <= Ls), 'atoms are outside the box');

    u = mean(vnew, 1);

    T_inst = m*sum(vecnorm(vnew, 2, 2).^2-norm(u)^2, 1)/(3*N);

    if t < STEPS_thermostat
       vnew = vnew*sqrt(Ts/T_inst); 
    end
