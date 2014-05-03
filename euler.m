function output = euler(TPitch, TYaw, TRoll)
    alpha = TPitch;
    beta = TYaw;
    gamma = TRoll;
    Rx = [1,0,0;                    0,cos(alpha),-sin(alpha);       0,sin(alpha),cos(alpha)];
    Ry = [cos(beta),0,sin(beta);    0,1,0;                          -sin(beta),0,cos(beta)];
    Rz = [cos(gamma),-sin(gamma),0; sin(gamma),cos(gamma),0;        0,0,1];
    output = Rx*Ry*Rz;
end

