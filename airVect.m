%outputs Air Vector as unit vector
function output = airVect(V_x_inertial, V_y_inertial, V_z_inertial, wind)

totalAir = [(wind(1) - V_x_inertial);
            (wind(2) - V_y_inertial);
            (wind(3) - V_z_inertial)];
V_x_inertial;
V_y_inertial;
V_z_inertial;
mag = sqrt(totalAir(1)^2 + totalAir(2)^2 + totalAir(3)^2);
if(mag==0)
    output = totalAir;
else
    output = totalAir/mag;
end
end