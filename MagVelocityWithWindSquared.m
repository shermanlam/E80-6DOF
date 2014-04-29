function ret=MagVelocitywithWindsquared(V_x_inertial, V_y_inertial, V_z_inertial, wind)
    ret = (wind(1) - V_x_inertial)^2 + (wind(2) - V_y_inertial)^2 + (wind(3) - V_z_inertial)^2;
end