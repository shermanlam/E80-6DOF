%Developing a 6DOF model for our rocket, accounting for inertial forces on
%the rocket, eventually accounting for gyroscoping forces.
%Our rocket is the Arreaux by Aerotech
clear all

%~~~~~***CONTROLS***~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TimeEnd = 10;               %End time of rocket flight in seconds
numPlotPoints = 800*TimeEnd;          %Number of points used to plot the final trajectory


%~~~~~~CONSTANTS & quantities~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%------Of Gyro
Omega_gyro = 43000/60*(2*pi); %Now in Rad/sec
mgyro = .300;
rgyro = .017145;
Lgyro = .5*mgyro*rgyro^2;

%~~~~~~Of Rocket
m = 1.422;          %Mass of with motor
L = 1.4;                    %Length of rocket~ FIRST ERROR
diam = 0.0483;
Area = 0.06762;              %Planform Area of the rocket
Area2 = 0.00183224752;
CG = 0.796;                 %Dist. from tip to center of pressure
CP = 1.14;                 %Dist. from tip to center of Gravity
CptoCG = CP -CG;            %radius distance between centre of press and cg

%~~~~~~Of the World (assuming constant)
DensityAirElevFix = 1.1839;
cc = 1.3;                   %Damping coefficient of rotation N/rad s
g = 9.81;                   %Gravitational Constant
uu = 1.814*10^(-5);         %Viscosity of Air N s/m^2 is constant
wind = [1;1;0];             %Direction and speed of wind, m/s inertial frame

LRollZ = 4.335*(10^(-4));     %These are the moments of inertia values from OpenRocket
LPitch = 0.21;
LYaw = 0.21;

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
                870;
                0;
                0;
                0;
                0; 
                0;
                0;
                0;
                0;
                0];
% initialConds = [-39.2215;
%                 6.3393;
%                 3502.8;
%                 -.4049;
%                 7.901;
%                 -2.7122;
%                 -.2344;
%                 -.1856;
%                 -1322.3;
%                 -1.8367;
%                 .5682;
%                 59.4505];
            

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

%     plott = linspace(0,6,100);
%     plot(plott, CL(plott))
    %kosher.
    
%~~~~~~FORCES IN BODY FRAME~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%drag force intermediate helper equations 
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
        abs(0.5*CL(AoA(V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, wind))...
        *DensityAirElev(Z_inertial,0)*Area*cos(AoA(V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, wind))*MagVelocityWithWindSquared(V_x_inertial, V_y_inertial,V_z_inertial,wind))...
        .*LiftNormalizeDirection(V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, wind);
                            

FdBody =@(Z_inertial,V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw) ...
        0.5*DensityAirElev(Z_inertial,0)*(Area2+Area*sin(AoA(V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, wind)))*...
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
        FdBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw) +...
        FThrustBody(t);
    
%~~~~~NEWTONS EQNS OF MOTION INTO INERTIAL FRAME~~~~~~~~~~~~~~~~~~~~~~~~~~~
LinAccInertial =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, TRoll, WPitch, WYaw, WRoll, t) ...
                  ((euler(-TPitch, -TYaw, -TRoll)\SumForceBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, TRoll, t)...
                    + m*cross([WPitch; WYaw; WRoll], [V_x_inertial; V_y_inertial; V_z_inertial]))...
                  +[0; 0; -m*g]+FLiftBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw))./m;

% LinAccInertial =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, TRoll, WPitch, WYaw, WRoll, t) ...
%                   (euler2(SumForceBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, TRoll, t),TPitch, TYaw, TRoll)...
%                     + m*cross([WPitch; WYaw; WRoll], [V_x_inertial; V_y_inertial; V_z_inertial])...
%                   +[0; 0; -m*g]+FLiftBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw))./m;
    
%~~~~~TORQUES IN THE BODY FRAME~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%normalized equations to determine direction of torque from gyro precession 
%in the Body Frame ARE NOW INSIDE THEIR OWN SEPARATE FUNCTIONS
%dPitch
%dYaw

%Now all the torque equations in the body frame

TorqueLift = @(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw)...
              CptoCG*cos(abs(AoA(V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, wind))) *cross(rocketVect(TPitch, TYaw), FLiftBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw));
            
TorqueGyro =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, WPitch, WYaw, t)...
              -0*(t<=8)*Lgyro*Omega_gyro*[WPitch; WYaw; 0];%-TorqueLift(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw)...
              %/norm(TorqueLift(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw));   %Do we have direction on this correct?

TorqueGyroRoll = 0*Omega_gyro * Lgyro/LRollZ * [0;0;-1];

TorqueDamping = @(WPitch, WYaw)...
               [-PitchDamping(WPitch)*sign(WPitch);-PitchDamping(WYaw)*sign(WYaw);0];   %Torque Damping only happens in Pitch angle direction.

%Putting all the T  werks together.
SumTorquesBody =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial,TPitch, TYaw, WPitch, WYaw, t)...
                TorqueDamping(WPitch, WYaw) + TorqueGyroRoll + TorqueGyro(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, WPitch, WYaw, t);
                %Gyro torque is added because it already accounts for
                %opposing the direction of motion, and is thus already
                %negative

              
%~~~~~ANGULAR MOMENTUMS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


ITensor = [LPitch, 0, 0;    %The tensor matrix. ooooh.
           0, LYaw, 0;
           0, 0, LRollZ];

%InverseTensor = inv(ITensor);  %Adding the inverse of the tensor.      

AngAccInertial =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, TRoll, WPitch, WYaw, WRoll, t)...
    ((euler(-TPitch, -TYaw, -TRoll)\SumTorquesBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial,TPitch, TYaw, WPitch, WYaw, t)...
    + cross([WPitch; WYaw; WRoll],[LPitch;LYaw;LRollZ].*[WPitch; WYaw; WRoll])...
    + TorqueLift(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw)) ./ [LPitch;LYaw;LRollZ]).*[1;1;0];
% 
% AngAccInertial =@(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw, TRoll, WPitch, WYaw, WRoll, t)...
%     ITensor\((euler2(SumTorquesBody(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial,TPitch, TYaw, WPitch, WYaw, t), TPitch, TYaw, TRoll))+...
%     cross([WPitch; WYaw; WRoll],ITensor*[WPitch; WYaw; WRoll])...
%     + TorqueLift(Z_inertial, V_x_inertial, V_y_inertial, V_z_inertial, TPitch, TYaw)) .* [1;1;0];

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
                indexat(LinAccInertial(r(3),r(4),r(5),r(6),r(7),r(8), r(9), r(10),r(11),r(12),t),1);
                indexat(LinAccInertial(r(3),r(4),r(5),r(6),r(7),r(8), r(9), r(10),r(11),r(12),t),2);
                indexat(LinAccInertial(r(3),r(4),r(5),r(6),r(7),r(8), r(9), r(10),r(11),r(12),t),3);
                r(10);
                r(11);
                r(12);
                indexat(AngAccInertial(r(3),r(4),r(5),r(6),r(7),r(8),r(9),r(10),r(11),r(12), t),1);
                indexat(AngAccInertial(r(3),r(4),r(5),r(6),r(7),r(8),r(9),r(10),r(11),r(12), t),2);
                indexat(AngAccInertial(r(3),r(4),r(5),r(6),r(7),r(8),r(9),r(10),r(11),r(12), t),3)];
 
%Using ODE45
%options = odeset('RelTol', .000001, 'AbsTol', .000001, 'Stats', 'on'); 
Solution= ode45(diffeq,[0,TimeEnd],initialConds);%, options);
% Solution= ode5(diffeq,linspace(0,TimeEnd,numPlotPoints),initialConds);

%PLOTTING
figure(1)
clf
plott = linspace(0,TimeEnd,numPlotPoints);         % time points for plotting

%Making points for all dimensions in the inertial frame.
Xspace = deval(Solution, plott, 1);
Yspace = deval(Solution, plott, 2);
Zspace = deval(Solution, plott, 3);
% Xspace = Solution(:,1);
% Yspace = Solution(:,2);
% Zspace = Solution(:,3);

%Ploting individual points, varying colour
 plot3(real(Xspace), real(Yspace), real(Zspace), 'Color', [1 0 0], 'Marker', 'o')
