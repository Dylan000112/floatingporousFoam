%% Moody model file for 4 catenary chains

% script for case setup %
dimensionNumber = 3;
waterLevel = 0.65;         % [m]       z-coordinate of mean water level
waterDensity = 1000.0;     % [kg/m??]   Density of water
airDensity = 1.0;         % [kg/m??]   Density of air
gravity = 1;

time.start  = 0;
time.end = 40;
time.cfl = 0.9;
time.scheme = 'RK3';

print.dt = 0.1; % output interval.

%% extra quadpoints used for increased ground contact performance.
% it makes the difference between stable and non-stable results. 
numLib.qPointsAdded = 10;

%-------------------------------------------------------------------------%
%--------------------------- Ground model input  -------------------------%
%-------------------------------------------------------------------------%
ground.type = 'springDampGround';
ground.level = -0.5;
ground.dampingCoeff = 1.0;
ground.frictionCoeff = 0.1;
ground.vc = 0.01;
ground.stiffness = 300.0e6;

%-------------------------------------------------------------------------%
%---------------------------- Type definition ----------------------------%
%-------------------------------------------------------------------------%
cableType1.diameter = 0.02;
cableType1.gamma0 = 1.00;
cableType1.CDn = 1.6;
cableType1.CDt = 0.05;
cableType1.CMn = 2.0;
cableType1.materialModel.type = 'biLinear';
cableType1.materialModel.EA =  3E4;
        
%-------------------------------------------------------------------------%
%------------------------------- Geometry --------------------------------%

vertexLocations = {
                       1    [-0.27      0.4      0.0 ];
                       2    [-0.25      0.4      0.43];
                       
                       3    [-0.27     -0.4      0.0 ];
                       4    [-0.25     -0.4      0.43];
                       
                       5    [ 0.27      0.4      0.0 ];
                       6    [ 0.25      0.4      0.43];
                       
                       7    [ 0.27     -0.4      0.0  ];
                       8    [ 0.25     -0.4      0.43];
                       
                            
                  };
                   
%                 9    [ 0.0      0.0      -0.016]; 

                         
%----- Object definitions -----%
cable1.typeNumber = 1;
cable1.startVertex = 1; %
cable1.endVertex = 2; %
%cable1.length = 0.43; % 
cable1.IC.type = "PreStrain";%'PreStrain';
cable1.IC.eps0 = 1.2E-3;
cable1.N =2; % 
% 
% Copy remaining info from cable1. short hand
cable2=cable1;
cable2.startVertex = 3; 
cable2.endVertex = 4; 
%
cable3=cable1;
cable3.startVertex = 5; 
cable3.endVertex = 6; 
%
cable4=cable1;
cable4.startVertex = 7; 
cable4.endVertex = 8; 

%-------------------------------------------------------------------------%
%                           Boundary conditions                           %
%-------------------------------------------------------------------------%

% Four anchors defined by vertexLocations
bc1.vertexNumber = 1;
bc1.type = 'dirichlet';
bc1.mode = 'fixed';
%
bc2=bc1;
bc2.vertexNumber = 3;
% 
bc3=bc1;
bc3.vertexNumber = 5;
% 
bc4=bc1;
bc4.vertexNumber = 7;

bc5.vertexNumber = 2;
bc5.type = 'dirichlet';
bc5.mode = 'externalPoint';

bc6=bc5;
bc6.vertexNumber = 4;

bc7=bc5;
bc7.vertexNumber = 6;

bc8=bc5;
bc8.vertexNumber = 8;

% --- API conectivitity --- %
API.bcNames = {'bc5','bc6','bc7','bc8'};
API.reboot= 'no';
API.syncOutput = 1; 
API.staggerTimeFraction= 0.5;


% To remove initial fluctuations due to very stiff ground
% staticSolver.relax = 1;
% staticSolver.relaxFactor = 0.9;

% ===== END OF FILE ===== %
