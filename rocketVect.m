%outputs rocket direction of pointing in unit vector form
function output = rocketVect(TPitch, TYaw)

output = [cos(TYaw)*sin(TPitch);
          sin(TYaw)*sin(TPitch);
          cos(TPitch)];
end