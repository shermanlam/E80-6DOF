function output = euler2(vect, TPitch, TYaw, TRoll )
    v1 = vect(1);
    v2 = vect(2);
    v3 = vect(3);
    out1 = v1*cos(TYaw)*cos(TRoll) + v2*cos(TPitch)*sin(-TRoll) + v3*sin(-TYaw)*cos(TPitch);
    out2 = v2*cos(TPitch)*cos(TRoll)+v1*sin(-TRoll)*cos(TYaw)+v3*sin(-TPitch)*cos(TYaw);
    out3 = v3*cos(TPitch)*cos(TYaw)+v1*sin(-TYaw)*cos(TRoll)+ v2*cos(TRoll)*sin(-TPitch);
    output = [out1;out2;out3];
end

