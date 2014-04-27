%density of air as function of altitude in meters above SEA LEVEL
%choose if you want a temperature. In CELCIUS. Otherwise leave zero.
%If left zero, will simply use atmospheric temperature model

%information from: http://en.wikipedia.org/wiki/Density_of_air
%based on standard atmospheric model
function output = DensityAirElev(height,myTemperature)

p0 = 101.325;               %sea level standard atmos. press kPA
T0 = 288.15;                %sea level standard temperature
g = 9.80665;                %earth-surface gravitational acceleration m/s^2
L = 0.0065;                 %temperature lapse rate,  K/m
R = 8.31447;                %ideal (universal) gas constant, J/(mol·K)
M = 0.0289644;              %molar mass of dry air kg/mol

%only include this temperature if you want 
if myTemperature <= 0
    Temperature = T0 - L*height;
else
    'Using your Own Temperature'
    Temperature = myTemperature + 273;
end

Pressure = 1000*(p0*(1 - (L*height)/T0 )^((g*M)/(R*L)));

%and so the density can be found from the ideal gas law as:
output = (Pressure*M)/(R*Temperature);