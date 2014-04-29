%This is our big angle of attack normalization function.
function output = AoA(V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, wind)

Vmag  = sqrt((wind(1) - V_x_inertial)^2 + (wind(2) - V_y_inertial)^2 + (wind(3) - V_z_inertial)^2);

if(Vmag == 0)
         output = 0;
else
    Runit = [sin(TPitch)*cos(TYaw);sin(TPitch)*sin(TYaw);cos(TPitch)];
    Vair  = [-(wind(1) - V_x_inertial); -(wind(2) - V_y_inertial); -(wind(3) - V_z_inertial)];
    Vair  = Vair/Vmag;
    output = abs(acos(dot(Runit, Vair)));
end
% This comment is pointless
end
