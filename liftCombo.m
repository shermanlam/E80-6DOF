function ret=liftCombo(V_x_inertial, V_y_inertial, V_z_inertial, wind)
    if((V_x_inertial^2 + V_y_inertial^2) == 0)
    output = [0;0;0];
    else
    output = [-V_x_inertial/sqrt(V_x_inertial^2 + V_y_inertial^2);...
              -V_y_inertial/sqrt(V_x_inertial^2 + V_y_inertial^2);...
              0];
    end
    ret =  output * ((wind(1) - V_x_inertial)^2 + (wind(2) - V_y_inertial)^2 + (wind(3) - V_z_inertial)^2); 
end