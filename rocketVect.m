%outputs rocket direction of pointing in unit vector form
function output = rocketVect(TPitch, TYaw, TRoll)
mag = sqrt(TPitch^2 + TYaw^2 + TRoll^2);

output = [TPitch;
          TYaw;
          TRoll]/mag;
end