function z_pred = ekf_h(x, params)
    % x      : [phi; theta; psi] (Roll, Pitch, Yaw)
    % params : struct containing g (gravity) and m_ref (inertial magnetic vector)
    
    phi   = x(1);
    theta = x(2);
    psi   = x(3);
    
    % Reference vectors in Inertial Frame (NED)
    % Standard gravity is typically [0; 0; 1] (normalized)
    % Magnetic reference is typically [mx; 0; mz] (normalized, accounts for dip angle)
    g_I = [0; 0; 9.81];
    
    % Pre-calculate trig terms
    sr = sin(phi);   cr = cos(phi);
    sp = sin(theta); cp = cos(theta);
    sy = sin(psi);   cy = cos(psi);

    % Rotation Matrix: Inertial to Body (R_x * R_y * R_z)
    R_IB = [ cp*cy,                cp*sy,               -sp;
             sr*sp*cy - cr*sy,     sr*sp*sy + cr*cy,     sr*cp;
             cr*sp*cy + sr*sy,     cr*sp*sy - sr*cy,     cr*cp ];

    % Predicted measurements
    h_acc = R_IB * g_I;
    h_mag = R_IB * params;

    % Combine into a single 6x1 measurement vector
    z_pred = [h_acc; h_mag];
end