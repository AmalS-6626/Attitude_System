function F = ekf_F_jac(x, u)
    % x = [phi; theta; psi]
    % u = [p; q; r; dt; delta]
    
    phi   = x(1);
    theta = x(2);
    p = u(1); q = u(2); r = u(3);
    dt = u(4);
    delta = u(5);

    % Handle singularity as in state propagation
    % if abs(theta) >= pi/2
    %     theta = sign(theta) * (pi/2 - delta);
    % end

    % Pre-calculate trig terms for efficiency
    sr = sin(phi);   cr = cos(phi);
    sp = sin(theta); cp = cos(theta);
    tp = sp/cp;      % tan(theta)
    cp2 = cp^2;      % cos^2(theta)

    % Continuous-time Jacobian (df/dx)
    dfdx = zeros(3,3);

    % Row 1: d(phi_dot)/dx
    dfdx(1,1) = q*cr*tp - r*sr*tp;          % d/d_phi
    dfdx(1,2) = (q*sr + r*cr) / cp2;        % d/d_theta
    dfdx(1,3) = 0;                          % d/d_psi

    % Row 2: d(theta_dot)/dx
    dfdx(2,1) = -q*sr - r*cr;               % d/d_phi
    dfdx(2,2) = 0;                          % d/d_theta
    dfdx(2,3) = 0;                          % d/d_psi

    % Row 3: d(psi_dot)/dx
    dfdx(3,1) = (q*cr - r*sr) / cp;         % d/d_phi
    dfdx(3,2) = (q*sr + r*cr) * sp / cp2;   % d/d_theta
    dfdx(3,3) = 0;                          % d/d_psi
 
    % Discrete-time Jacobian: F = I + df/dx * dt
    F = eye(3) + dfdx * dt;
end