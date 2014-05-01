%Developing a 6DOF model for our rocket, accounting for inertial forces on
%the rocket, eventually accounting for gyroscoping forces.
%Our rocket is the Arreaux by Aerotech
clear all

%~~~~~***CONTROLS***~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TimeEnd = 5;               %End time of rocket flight in seconds
numPtsPlot = 100*TimeEnd;          %Number of points used to plot the final trajectory


%~~~~~~CONSTANTS & quantities~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%------Of Gyro
Omega_gyro = 1;
mgyro = .500;
rgyro = .02286;
Lgyro = .5*mgyro*rgyro^2;

%~~~~~~Of Rocket
m = 0.449 + mgyro;          %Mass of with motor
L = 1.4;                    %Length of rocket~ FIRST ERROR
diam = 0.0483;
Area = 0.0073;              %Planform Area of the rocket
CG = 0.746;                 %Dist. from tip to center of pressure
CP = 0.942;                 %Dist. from tip to center of Gravity
CptoCG = CP -CG;            %radius distance between centre of press and cg

%~~~~~~Of the World (assuming constant)
DensityAirElevFix = 1.1839;
cc = 1.3;                   %Damping coefficient of rotation N/rad s
g = 9.81;                   %Gravitational Constant
uu = 1.814*10^(-5);         %Viscosity of Air N s/m^2 is constant
wind = [1;0;0];             %Direction and speed of wind, m/s inertial frame


%~~~~~INITIAL CONDITIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Where
% r(1) = x                      
% r(2) = y
% r(3) = z                      %Above Sea Level
% r(4) = xdot
% r(5) = ydot
% r(6) = zdot
% r(7) = thetax Tpitch
% r(8) = thetay Tyaw
% r(9) = thetaz Troll
% r(10) = thetadotx WPitch
% r(11) = thetadoty WYaw
% r(12) = thetadotz Wroll
initialConds = [0;
                0;
                3000;
                0.1;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0];
            

%~~~~~~IMPORTING FOR SPLINE CURVES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Thrust Curve Importing:
    thrustfile = 'G75.csv';     %Specify CSV file where 1st col = time, 
                                %and 2nd col = thrust in newtons
    thrustRaw = csvread(thrustfile);
    ThrustTimes = thrustRaw(:,1);
    ThrustForceN = thrustRaw(:,2);
%Coefficient of Drag Importing, as function of reynolds number
    dragfile = 'G75cd2.csv';    %Specify CSV file where 1st col = Reynolds Num, 
                                %and 2nd col = Coefficient of Drag
    dragRaw = csvread(dragfile);
    ReynoldsCD = dragRaw(:,1);
    CDsvsReynold = dragRaw(:,2);
%CL vs Angle of Attack Curve Importing:
    CLfile = 'G75AoAvsCL3.csv';     %Specify CSV file where 1st col = ANGLE, 
                                    %and 2nd col = CL which is unitless
    CLRaw = csvread(CLfile);
    CLAngles = CLRaw(:,1);
    CLVals = CLRaw(:,2);
%RotationalDamping vs WPitch (Rad/sec) Curve Importing:
    PitchDampFile = 'G75pitchdamp.csv';     %Specify CSV file where 1st col = Damping, 
                                            %and 2nd col = WPitch, in
                                            %Rad/second
    PitchDampRaw = csvread(PitchDampFile);
    PitchRotRates = PitchDampRaw(:,1);
    DampingPitchVal = PitchDampRaw(:,2);
    

%~~~~~CREATING SPLINE CURVES and their EQUATIONS ~~~~~~~~~~~~~
%Creating a Spline Curves for the various more complicated parameters
    ThrustSpline = spline(ThrustTimes, ThrustForceN);
    ThrustCurve =@(t) ppval(ThrustSpline, t);
%Now for coefficient of drag and equation
    CdSpline = spline(ReynoldsCD, CDsvsReynold);
    Cd =@(Re) ppval(CdSpline, Re);
%Now for coefficient of LIFT and its equation vs. total magnitude of Angle
%of Attack against true airspeed and its function
    CLSpline = spline(CLAngles, CLVals);
    CL =@(ThetaMag) ppval(CLSpline, ThetaMag);
%Now for Spline coefficient of Pitch damping of the rocket, and creating an
%equation for Damping vs. WPitch
    PitchDampSpline = spline(PitchRotRates, DampingPitchVal);
    PitchDamping =@(WPitch) ppval(PitchDampSpline, WPitch);

    %kosher.
    
%~~~~~~FORCES IN BODY FRAME~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%drag force intermediate helper equations 
Magnitude =@(x,y,z) sqrt(x^2 + z^2 + y^2);
Reynolds =@(V_air,Z_inertial) (DensityAirElev(Z_inertial,0)*V_air*L)/uu;
                        
%AoA moved to separate function. Must include wind now.

%these first two functions have the incident wind accounted for,
%All the force equations in body frame:
% FLiftBody =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw) ...
%              0.5*CL(AoA(V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, wind))...
%              *DensityAirElev(Z_inertial,0)*Area*...
%              MagVelocitywithWindsquared(V_x_inertial, V_y_inertial, V_z_inertial)...
%              *LiftNormalizeDirection(V_x_inertial, V_y_inertial); 

% FLiftBody =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw) ...
%              0.5*CL(AoA(V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, wind))...
%              *DensityAirElev(Z_inertial,0)*Area*...
%              liftCombo(V_x_inertial, V_y_inertial, V_z_inertial, wind); 
FLiftBody =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw)...
        0.5*CL(AoA(V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, wind))...
        *DensityAirElev(Z_inertial,0)*Area*MagVelocityWithWindSquared(V_x_inertial, V_y_inertial,V_z_inertial,wind)...
        .*LiftNormalizeDirection(V_x_inertial, V_y_inertial, V_z_inertial, wind);
                            

FdBody =@(Z_inertial,V_x_inertial, V_y_inertial, V_z_inertial, TPitch) ...
        0.5*DensityAirElev(Z_inertial,0)*Area*...
        MagVelocityWithWindSquared(V_x_inertial, V_y_inertial, V_z_inertial,wind)*cos(TPitch)^2*...
        Cd(Reynolds(sqrt(MagVelocityWithWindSquared(V_x_inertial,V_z_inertial, V_y_inertial, wind)), Z_inertial))...
        .*[0;0;-1];    %Drag

% FdBody =@(Z_inertial,V_x_inertial, V_y_inertial, V_z_inertial, TPitch) ...
%         0.5*DensityAirElevFix*Area*...
%         MagVelocitywithWindsquared(V_x_inertial, V_y_inertial, V_z_inertial)*cos(TPitch)^2*...
%         Cd(Reynolds(sqrt(MagVelocityWithWindsquared(V_x_inertial,V_z_inertial, V_y_inertial,wind)), Z_inertial))...
%         .*[0;0;-1];    %Drag

% FdBody =@(Z_inertial,V_x_inertial, V_y_inertial, V_z_inertial, TPitch) [0;0;0];
% FLiftBody =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw) [0;0;0];
    
FThrustBody =@(t) [0; 0; ThrustCurve(t)];

%Putting all dem forces together, in the Body Frame!
SumForceBody =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, TRoll, t)...
        FdBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch) +...
        FThrustBody(t);
    
%~~~~~NEWTONS EQNS OF MOTION INTO INERTIAL FRAME~~~~~~~~~~~~~~~~~~~~~~~~~~~
LinAccInertial =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, TRoll, WPitch, WYaw, WRoll, t) ...
                  ((SumForceBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, TRoll, t)...
                    + m*cross([WPitch; WYaw; WRoll], [V_x_inertial; V_y_inertial; V_z_inertial]))...
                  +[0; 0; -m*g]+FLiftBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw))./m;

% LinAccInertial =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, WPitch, WYaw, WRoll, t) ...
%                   (SumForceBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, t)...
%                   +[0; 0; -m*g])./m;  
    
%~~~~~TORQUES IN THE BODY FRAME~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%normalized equations to determine direction of torque from gyro precession 
%in the Body Frame ARE NOW INSIDE THEIR OWN SEPARATE FUNCTIONS
%dPitch
%dYaw

%Now all the torque equations in the body frame
TorqueGyro =@(WPitch, WYaw)...
           (Omega_gyro * Lgyro).*[dPitch(WPitch, WYaw);dYaw(WPitch, WYaw);0];   %Do we have direction on this correct?

TorqueLift = @(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw)...
            CptoCG.*FLiftBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw);

TorqueDamping = @(WPitch)...
                PitchDamping(WPitch)*WPitch.*[1;0;0];   %Torque Damping only happens in Pitch angle direction.

%Putting all the T  werks together.
SumTorquesBody =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial,TPitch, TYaw, WPitch, WYaw)...
                TorqueDamping(WPitch) + TorqueGyro(WPitch, WYaw)... 
                - TorqueLift(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw);
                %Gyro torque is added because it already accounts for
                %opposing the direction of motion, and is thus already
                %negative


              
%~~~~~ANGULAR MOMENTUMS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LRollZ = 2.9*(10^(-4));     %These are the moments of inertia values from OpenRocket
LPitch = 0.12994;
LYaw = 0.12994;

ITensor = [LPitch, 0, 0;    %The tensor matrix. ooooh.
           0, LYaw, 0;
           0, 0, LRollZ];

AngAccInertial =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, TRoll, WPitch, WYaw, WRoll)...
    SumTorquesBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial,TPitch, TYaw, TRoll, WPitch, WYaw)+...
    cross([WPitch; WYaw; WRoll], (ITensor*[TPitch; TYaw; TRoll]));



%~~~~~~INDEXAT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
indexat =  @(vector,indices) vector(indices);


%~~~~SOLVING IT WITH ODE45~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% r(1) = x
% r(2) = y
% r(3) = z
% r(4) = xdot
% r(5) = ydot
% r(6) = zdot
% r(7) = thetax Tpitch
% r(8) = thetay Tyaw
% r(9) = thetaz Troll
% r(10) = thetadotx WPitch
% r(11) = thetadoty WYaw
% r(12) = thetadotz Wroll
diffeq =@(t,r) [r(4);
                r(5);
                r(6);
                indexat(LinAccInertial(r(3),r(4),r(5),r(6),r(7),r(8),r(10),r(11),r(12),t),1);
                indexat(LinAccInertial(r(3),r(4),r(5),r(6),r(7),r(8),r(10),r(11),r(12),t),2);
                indexat(LinAccInertial(r(3),r(4),r(5),r(6),r(7),r(8),r(10),r(11),r(12),t),3);
                r(10);
                r(11);
                r(12);
                indexat(AngAccInertial(r(3),r(4),r(5),r(6),r(7),r(8),r(9),r(10),r(11),r(12)),1);
                indexat(AngAccInertial(r(3),r(4),r(5),r(6),r(7),r(8),r(9),r(10),r(11),r(12)),2);
                indexat(AngAccInertial(r(3),r(4),r(5),r(6),r(7),r(8),r(9),r(10),r(11),r(12)),3)];
 
%Using ODE45
%options = odeset('RelTol', .000001, 'AbsTol', .000001, 'Stats', 'on'); 
%Solution= ode45(diffeq,[0,TimeEnd],initialConds);%, options);
Solution= ode5(diffeq,linspace(0,TimeEnd,1000),initialConds);

%PLOTTING
figure(1)
clf
plott = linspace(0,TimeEnd,1000);         % time points for plotting

%Making points for all dimensions in the inertial frame.
% Xspace = deval(Solution, plott, 1);
% Yspace = deval(Solution, plott, 2);
% Zspace = deval(Solution, plott, 3);
Xspace = Solution(:,1);
Yspace = Solution(:,2);
Zspace = Solution(:,3);

%Ploting individual points, varying colour
 plot3(Xspace, Yspace, Zspace, 'Color', [1 0 0], 'Marker', 'o')
