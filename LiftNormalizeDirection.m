%Correction function for normaization of lift direction.
%Last matrix with negative signs accounts for opposing of direction of rotation

function output = LiftNormalizeDirection(V_x_inertial, V_y_inertial)

if((V_x_inertial^2 + V_y_inertial^2) == 0)
    output = [0;0;0];
else
    output = [-V_x_inertial/sqrt(V_x_inertial^2 + V_y_inertial^2);...
              -V_y_inertial/sqrt(V_x_inertial^2 + V_y_inertial^2);...
              0];
end
end