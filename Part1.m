%part1.m - version modified of DemoEKF.m originally provided by lecturer,
%MTRN 4010 - 2014 by Karan Narula z3435335
%Project 1 part 1
function main()
stdDevGyro = 2*pi/180 ;          % 2 degrees/second , standard deviation of the gyros' noise
stdDevSpeed = 0.3 ;   % speed sensor standard deviation = 0.3m/sec 
% ... errors in the range measurements (25cm, standard.dev.)
sdev_rangeMeasurement = 0.25 ;          % std. of noise in range measurements. 0.25m
sdev_bearingMeasurement = deg2rad(5);  % std of noise in the bearing measurements. 5 degrees (high value !)
%the actual laser scanner should not flactuate more than 1-2 degrees

Dt=0.05 ;                       % "sample time", 50ms
Li = 5000 ;                     % number of iterations of experiment 
DtObservations=0.250 ;          % laser "sample time" of observations, 4Hz, approximately

% Number of landmarks in use  
n_usedLanmarks = 3 ;    
% just to create some landmarks
global NavigationMap;
NavigationMap = CreateSomeMap(n_usedLanmarks) ;  %creates a artificial map!
%initial belief about our state variables
Xe = [ 0; 0;pi/2 ] ;        
P = zeros(3,3) ;            % absolute confidence in initial belief

% Some buffers to store the intermediate values during the experiment
Xreal_History= zeros(3,Li) ;
Xe_History= zeros(3,Li);
%uncertainty in the process model x-> 1 cm^2; y-> 1cm^2; heading -> 1 degree
Q1 = diag( [ (0.01)^2 ,(0.01)^2 , (1*pi/180)^2]) ;

%co-variance of the inputs, a diagonal matrix since the noise of gyro and speed is independent
Pu = diag([stdDevSpeed*stdDevSpeed, stdDevGyro*stdDevGyro]);


time=0 ;
% init simulator of process.
InitSimulation(stdDevSpeed,stdDevGyro,sdev_rangeMeasurement,sdev_bearingMeasurement, DtObservations);
for i=1:Li,     % loop
        
    time = time+Dt ;    
    SimuPlatform(time,Dt);     
    %get the input that has been simulated and polluted with noise
    [Noisy_speed,Noisy_GyroZ]=GetProcessModelInputs();
    %prediction using process model. Expected value of Xe @ k+1 given k
    Xe    = RunProcessModel(Xe,Noisy_speed,Noisy_GyroZ,Dt) ;
    
    % Estimate new covariance after prediction
    J = [ [1,0,-Dt*Noisy_speed*sin(Xe(3))  ]  ; [0,1,Dt*Noisy_speed*cos(Xe(3))] ;    [ 0,0,1 ] ] ;
    %3x2 jacobian matrix of process model wrt to input
    Ju = [Dt*cos(Xe(3)), 0; Dt*sin(Xe(3)), 0; 0, Dt];
    %uncertainty in the inputs
    Qu = Ju*Pu*Ju';
    %Total uncertainty contribution from both the inputs and the process Model
    Q = Q1 + Qu;    
    % then I calculate the new coveraince, after the prediction P(K+1|K) = J*P(K|K)*J'+Q ;
    P = J*P*J'+Q ;


    % .. Get range measuremens, if those are available.
    [nDetectedLandmarks,MasuredRanges,MeasuredBearing, IDs]=GetObservationMeasurements();
    %In this part, measured ranges will not be used
    
    if nDetectedLandmarks>0,     
        
       for u=1:nDetectedLandmarks
           ID = IDs(u);            % landmark ID
           %find expected bearing using information from prediction and current measurement
           theta = atan2(NavigationMap.landmarks(ID,2) - Xe(2),NavigationMap.landmarks(ID,1) - Xe(1));
           ExpectedBearing = theta - Xe(3) + pi/2;
           
           %here, we refer to lecture note to see the formula of H matrix
           %for bearing measurement
           eDX = (NavigationMap.landmarks(ID,1)-Xe(1)) ;      % (xu-x)
           eDY = (NavigationMap.landmarks(ID,2)-Xe(2)) ;      % (yu-y)
           eDD = eDX.*eDX + eDY.*eDY;
           H = [eDY/eDD, -eDX/eDD, -1];
           % Evaluate residual (innovation)  "Y-h(Xe)"
           z  = MeasuredBearing(u) - ExpectedBearing ;
           %clamp the value to be between -pi and pi
           z = mod(z+pi, 2*pi)-pi;
           %uncertainty in the measurement values
           R = sdev_bearingMeasurement*sdev_bearingMeasurement*4;
           S = R + H*P*H' ;
           iS = inv(S) ;               
           
           K = P*H'/S ;           % Kalman gain
           
           Xe = Xe+K*z ;           % update the  expected value
           P = P-P*(H'/S)*H*P ;     % update the Covariance
                 
       end      
     
    end;   
     Xreal_History(:,i) = GetCurrentSimulatedState() ;
     Xe_History(:,i)    = Xe ;
end ;    
      


SomePlots(Xreal_History,Xe_History,NavigationMap) ;

return ;          

    
function Xnext=RunProcessModel(X,speed,GyroZ,dt) 
    Xnext = X + dt*[ speed*cos(X(3)) ;  speed*sin(X(3)) ; GyroZ ] ;
return ;




function [ranges, Bearing, IDs] = GetMeasurementsFomNearbyLandmarks(X,map)
    dx= map.landmarks(:,1) - X(1) ;
    dy= map.landmarks(:,2) - X(2) ;
    ranges = sqrt((dx.*dx + dy.*dy)) ;
    %here, adding a statement to return the measured bearing 
    Bearing = atan2(dy, dx) - X(3) + pi/2;

    IDs = [1:map.nLandmarks];
return ;

function [speed,GyroZ] = SimuControl(X,t)
    speed = 2 ;                                         
    GyroZ = 3*pi/180 + sin(0.1*2*pi*t/50)*.02 ;         
return ;

function map = CreateSomeMap(n_used)
    landmarks = [  [ -40,0 ];[ -0,-20 ];[ 10,10 ] ;[ 30,10 ]  ] ;
    map.landmarks = landmarks(1:n_used,:) ;
    map.nLandmarks = n_used ;
return ;



function InitSimulation(stdDevSpeed,stdDevGyro,sdev_rangeMeasurement,sdev_bearingMeasurement, DtObservations)
    global ContextSimulation;
    ContextSimulation.Xreal = [ 0; 0;pi/2 ] ;     
    ContextSimulation.stdDevSpeed = stdDevSpeed;
    ContextSimulation.stdDevGyro = stdDevGyro;
    ContextSimulation.Xreal = [0;0;pi/2];
    ContextSimulation.speed=0;
    ContextSimulation.GyroZ=0;
    ContextSimulation.sdev_rangeMeasurement=sdev_rangeMeasurement; 
    ContextSimulation.sdev_bearingMeasurement = sdev_bearingMeasurement;
    ContextSimulation.DtObservations=DtObservations;
    ContextSimulation.timeForNextObservation= 0;
    ContextSimulation.CurrSimulatedTime=0;
return;


function [Noisy_speed,Noisy_GyroZ]=GetProcessModelInputs()
    global ContextSimulation;
    Noisy_speed =ContextSimulation.speed+ContextSimulation.stdDevSpeed*randn(1) ;
    Noisy_GyroZ =ContextSimulation.GyroZ+ContextSimulation.stdDevGyro*randn(1);
return;

function  SimuPlatform(time,Dt)
    global ContextSimulation;    
    [ContextSimulation.speed,ContextSimulation.GyroZ] = SimuControl(ContextSimulation.Xreal,time) ;    

    ContextSimulation.Xreal = RunProcessModel(ContextSimulation.Xreal,ContextSimulation.speed,ContextSimulation.GyroZ,Dt) ;
    ContextSimulation.CurrSimulatedTime = ContextSimulation.CurrSimulatedTime+Dt;
return;


function [nDetectedLandmarks,MasuredRanges, MeasuredBearing, IDs]=GetObservationMeasurements(map)
    global ContextSimulation NavigationMap;       
   
    if ContextSimulation.CurrSimulatedTime<ContextSimulation.timeForNextObservation,
        nDetectedLandmarks=0;
        MasuredRanges=[];
        MeasuredBearing = [];
        IDs=[];
        return ; 
    end;
        
    ContextSimulation.timeForNextObservation = ContextSimulation.timeForNextObservation+ContextSimulation.DtObservations;
    
        [RealRanges,RealBearing, IDs] = GetMeasurementsFomNearbyLandmarks(ContextSimulation.Xreal,NavigationMap) ;
              
        nDetectedLandmarks = length(RealRanges) ;
      
        noiseInMeasurements_range= ContextSimulation.sdev_rangeMeasurement*randn(size(RealRanges));
        noiseInMeasurements_bearing= ContextSimulation.sdev_bearingMeasurement*randn(size(RealBearing));
        
        MasuredRanges = RealRanges +  noiseInMeasurements_range ;
        MeasuredBearing = RealBearing + noiseInMeasurements_bearing;
 return;

 
    function X=GetCurrentSimulatedState()
        global ContextSimulation;
        X=ContextSimulation.Xreal;
        
        
     return;   

function SomePlots(Xreal_History,Xe_History,map) ;



figure(2) ; clf ; hold on ;
plot(Xreal_History(1,:),Xreal_History(2,:),'b') ;
plot(Xe_History(1,:),Xe_History(2,:),'r') ;
plot(map.landmarks(:,1),map.landmarks(:,2),'*r') ;
ii = [1:225:length(Xe_History(1,:))] ;
quiver(Xe_History(1,ii),Xe_History(2,ii),5*cos(Xe_History(3,ii)),5*sin(Xe_History(3,ii)),'r' ) ;
quiver(Xreal_History(1,ii),Xreal_History(2,ii),5*cos(Xreal_History(3,ii)),5*sin(Xreal_History(3,ii)),'b' ) ;
plot(Xe_History(1,ii),Xe_History(2,ii),'+r') ;
plot(Xreal_History(1,ii),Xreal_History(2,ii),'+b') ;
title('trajectories (blue:Real, red:EKF)') ;

figure(3) ; clf ; 
subplot(311) ; plot(Xreal_History(1,:)-Xe_History(1,:)) ;ylabel('x-xe (m)') ;
title('Performance EKF') ;
subplot(312) ; plot(Xreal_History(2,:)-Xe_History(2,:)) ;ylabel('y-ye (m)') ;
subplot(313) ; plot(180/pi*(Xreal_History(3,:)-Xe_History(3,:))) ;ylabel('heading error (deg)') ;

return ;