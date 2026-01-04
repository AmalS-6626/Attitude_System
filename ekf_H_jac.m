function H = ekf_H_jac(x, params)
    % x      : [phi; theta; psi]
    % params : [g; mx; my; mz]
    
    phi   = x(1);
    theta = x(2);
    psi   = x(3);
    
    g  = 9.81;
    mI = params(1:3);
    mN = mI(1); mE = mI(2); mD = mI(3);

    % Pre-calculate trig
    sr = sin(phi);   cr = cos(phi);
    sp = sin(theta); cp = cos(theta);
    sy = sin(psi);   cy = cos(psi);

    H = zeros(6,3);

    %% Row 1-3: Accelerometer Jacobians (dh_acc/dx)
    % dh_acc / d_phi
    H(1,1) = 0;
    H(2,1) = cr*cp*g;
    H(3,1) = -sr*cp*g;

    % dh_acc / d_theta
    H(1,2) =  -cp*g;
    H(2,2) =  -sr*sp*g;
    H(3,2) =  -cr*sp*g;

    % dh_acc / d_psi
    H(1,3) = 0;
    H(2,3) = 0;
    H(3,3) = 0;

    %% Row 4-6: Magnetometer Jacobians (dh_mag/dx)
    % dh_mag / d_phi
    H(4,1) = 0;
    H(5,1) = (cr*sp*cy + sr*sy)*mN + (cr*sp*sy - sr*cy)*mE + cr*cp*mD;
    H(6,1) = (-sr*sp*cy + cr*sy)*mN + (-sr*sp*sy - cr*cy)*mE - sr*cp*mD;

    % dh_mag / d_theta
    H(4,2) = -sp*cy*mN - sp*sy*mE - cp*mD;
    H(5,2) = sr*cp*cy*mN + sr*cp*sy*mE - sr*sp*mD;
    H(6,2) = cr*cp*cy*mN + cr*cp*sy*mE - cr*sp*mD;

    % dh_mag / d_psi
    H(4,3) = -cp*sy*mN + cp*cy*mE;
    H(5,3) = (-sr*sp*sy - cr*cy)*mN + (sr*sp*cy - cr*sy)*mE;
    H(6,3) = (-cr*sp*sy + sr*cy)*mN + (cr*sp*cy + sr*sy)*mE;
end