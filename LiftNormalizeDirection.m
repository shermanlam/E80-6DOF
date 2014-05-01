%Correction function for normaization of lift direction.
%Last matrix with negative signs accounts for opposing of direction of rotation

function output = LiftNormalizeDirection(V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, wind)
    airspeed2 = MagVelocityWithWindSquared(V_x_inertial, V_y_inertial, V_z_inertial, wind);
    if((sqrt(airspeed2) == 0))
        output = [0;0;0];
    else
        orientationVector = rocketVect(TPitch, TYaw);
        windVector = airVect(V_x_inertial, V_y_inertial, V_z_inertial, wind);
        output = windVector - dot(windVector, orientationVector).*orientationVector;
        if(norm(output) ~= 0)
            output = output ./ norm(output);
        end
    end
end