%% Supporting Functions
function [bondModel, springModel, crackModel] = getMaterialModels_TensionStiff(modelType)
%% GETMATERIALMODELS - Returns material model parameters based on analysis type
%
% INPUTS:
%   modelType - String specifying analysis type: 'RC', 'FRC', 'HRC', 'PAVEMENT'
%
% OUTPUTS:
%   bondModel - Structure containing bond-slip model parameters for rebar and fiber
%   springModel - Structure containing Friction model parameters/ Transverse yarn parameters / Transverse Reinforcement
%   crackModel - Structure containing Crack-Width model parameters for FRC, HRC, and PAVEMENT
%
% AUTHOR: Chidchanok Pleesudjai

switch upper(modelType)
    case 'RC'  % Reinforced Concrete
        [bondModel, springModel, crackModel] = getRC_Models();
        
    case 'FRC'  % Fiber Reinforced Concrete
        [bondModel, springModel, crackModel] = getFRC_Models();
        
    case 'HRC'  % Hybrid Reinforced Concrete
        [bondModel, springModel, crackModel] = getHRC_Models();
        
    case 'PAVEMENT'  % Pavement with friction
        [bondModel, springModel, crackModel] = getPavement_Models();
        
    otherwise
        error('Unknown model type. Use RC, FRC, HRC, or PAVEMENT');
end

end

function [bondModel, springModel, crackModel] = getRC_Models()
%% Standard Reinforced Concrete Model
% Bond-slip model for concrete-rebar interface
bondModel.slip = [0, 0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0] * 0.1;  % mm
bondModel.stress = [0, 8.5, 12.3, 14.2, 15.8, 12.4, 8.9, 5.2, 2.1];    % MPa
bondModel.failStress = 0;
bondModel.active = true;

% No spring model for standard RC
springModel.slip = [0, 1];
springModel.force = [0, 0];
springModel.failForce = 0;
springModel.active = false;

% No crack-width model for standard RC
crackModel.width = [0, 1];
crackModel.stress = [0, 0];
crackModel.active = false;
end

function [bondModel, springModel, crackModel] = getFRC_Models()
%% Fiber Reinforced Concrete Model  
% Fiber-matrix bond model
bondModel.slip = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0] * 0.1;
bondModel.stress = [0, 2.8, 4.2, 5.1, 4.8, 3.9, 2.7, 1.8, 0.9];
bondModel.failStress = 0;
bondModel.active = true;

% Crack bridging spring model (Transverse yarn parameters)
springModel.slip = [0, 0.05, 0.15, 0.5, 1.0] * 0.1;
springModel.force = [0, 15, 25, 18, 8];  % N/mm
springModel.failForce = 0;
springModel.active = false;

% Crack-width model for FRC (Stress-Crack width relationship)
crackModel.width = [0, 0.1, 0.3, 0.6, 1.0, 2.0];  % mm - Crack width
crackModel.stress = [0, 3.5, 2.8, 1.9, 1.2, 0.5]; % MPa - Residual tensile stress
crackModel.active = true;
end

function [bondModel, springModel, crackModel] = getHRC_Models()
%% Hybrid Reinforced Concrete Model (Rebar + Fibers)
% Enhanced bond model for hybrid system (Original parameters)
bondModel.slip = [0, 0.059331, 0.107162, 0.150086, 0.205674, 0.260493, ...
                  0.319289, 0.365804, 0.408920, 0.472656] * 0.1;
bondModel.stress = [0, 953.428509, 821.694853, 768.335547, 657.276218, ...
                    461.946413, 296.173598, 223.582923, 187.964676, 189.243689];
bondModel.failStress = 0;
bondModel.active = true;

% Enhanced spring model for transverse reinforcement (Original parameters)
springModel.slip = [0, 0.001, 0.15, 0.40, 1.00] * 0.0393701;
springModel.force = [0, 3, 2, 0.3, 0] * 0.01 * 145;  % Will be scaled by Ac in main program
springModel.failForce = 0;
springModel.active = false;

% Crack-width model for HRC (Enhanced due to hybrid reinforcement)
crackModel.width = [0, 0.05, 0.15, 0.4, 0.8, 1.5];  % mm - Crack width
crackModel.stress = [0, 4.2, 3.8, 2.5, 1.8, 1.0];   % MPa - Residual tensile stress
crackModel.active = true;
end

function [bondModel, springModel, crackModel] = getPavement_Models()
%% Pavement Model with Friction
% Modified bond model for pavement applications
bondModel.slip = [0, 0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0];
bondModel.stress = [0, 6.2, 8.9, 10.5, 9.8, 6.7, 3.2, 1.1];
bondModel.failStress = 0;
bondModel.active = true;

% Friction-slip model for pavement interface (Spring Model)
springModel.slip = [0, 0.02, 0.1, 0.3, 0.8, 2.0];  % mm
springModel.force = [0, 50, 80, 120, 100, 60];      % N/mm - Friction force
springModel.failForce = 10;
springModel.active = true;

% Crack-width model for pavement (Accounts for thermal and shrinkage cracking)
crackModel.width = [0, 0.2, 0.5, 1.0, 2.0, 4.0];   % mm - Crack width
crackModel.stress = [0, 2.1, 1.8, 1.2, 0.8, 0.3];  % MPa - Residual tensile stress
crackModel.active = true;
end

% Additional helper functions for the tension stiffening model would be
% implemented below this point, including:
% - Bond-slip relationship functions
% - Crack formation detection algorithms  
% - Stiffness matrix assembly routines
% - Solution convergence procedures