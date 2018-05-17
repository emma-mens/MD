function [rnew,vnew,T_inst]=take_one_step(N,r,v,F,t,dt,Ls,STEPS_thermostat,alpha,Ts)

    m = 1; 

    % calculate acceleration
    a = F/m;

    % update positions
    v = v + 0.5*a*dt;
    rnew = r + v*dt + 0.5*a*dt^2;

    % apply boundary condition on rnew
    rnew(rnew < 0 | rnew > Ls) = mod(rnew(rnew < 0 | rnew > Ls), Ls);

    u = mean(v, 1);

    T_inst = m*sum(vecnorm(v, 2, 2).^2-norm(u)^2, 1)/(3*N);

    if t < STEPS_thermostat
       v = v*sqrt(Ts/T_inst); 
    end
    
    vnew = v + 0.5*a*dt;
