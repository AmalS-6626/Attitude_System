function x_next = ekf_stateprop(x, u)
    phi = x(1); theta = x(2); psi = x(3);
    pqr = u(1:3); dt = u(4); delta = u(5);

    % Prevent singularity at pitch = +/- 90 degrees
    % if abs(theta) >= pi/2
    %     theta = sign(theta) * (pi/2 - delta);
    % end

    cp = cos(theta); 
    tp = tan(theta); % Simplified sp/cp
    sr = sin(phi);   cr = cos(phi);

    % Transformation matrix
    body_euler = [ 1, sr*tp, cr*tp; 
                   0, cr,    -sr; 
                   0, sr/cp, cr/cp ];

    x_dot = body_euler * pqr;
    x_next = x + x_dot * dt;
    
    % Wrap angles to [-pi, pi]
    % x_next = atan2(sin(x_next), cos(x_next)); 
end