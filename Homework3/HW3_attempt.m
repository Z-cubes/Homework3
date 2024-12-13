clear all
close all
clc


% Load data from Excel files
Years = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'UIF, Net Evap Data', 'Range', 'B4:B73'); % Years
I1 = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'UIF, Net Evap Data', 'Range', 'C4:C73'); % Annual Avg Flow (cfs)
E1_R = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'UIF, Net Evap Data', 'Range', 'C4:C73'); % Net Evap (inches/yr)
MAWS1 = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'DemandFunctions', 'Range', 'C6:C16'); % Marginal Value for City demand function ($/acre-ft)
MV1 = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'DemandFunctions', 'Range', 'D6:D16'); % Mean Annual Water Supply for City deman function (cfs)     
MAWS3 = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'DemandFunctions', 'Range', 'L4:L16'); % Marginal Value for Farmers demand function ($/acre-ft
MV3 = readtable('CEE4211 Homework 2 Data.xlsx', 'Sheet', 'DemandFunctions', 'Range', 'M4:M16'); % Mean Annual Water Supply for farmers deman function (cfs) 

% Convert tables to arrays  
Years = table2array(Years);
I1 = table2array(I1);
E1_R = table2array(E1_R);
MAWS1 = table2array(MAWS1);
MV1 = table2array(MV1);
MAWS3 = table2array(MAWS3);
MV3 = table2array(MV3);

% Initialize variable
U1original = 800; % determined target shares (cfs) Annual avg
U3original = 1100; % determined target shares (cfs) Annual avg
u1_R = U1original; % first iteration 
u3_R = U3original; % first iteration 
S1_R_original = 43.56; % Starting Storage (bcf)
S1_R = S1_R_original; % first iteration
S1_R_max = 43.56; % Maximum Storage (bcf) 
S1_R_min = 0; % Minimum Storage (bcf)
A1_R_original = 43.56e7; % Starting Reservoir Surface area (bsqf) for an assume staring depth of 100ft
A1_R = A1_R_original;  % assumed constant area 
Sp1_R_original = 0; % Staring spillage (bcf)
Sp1_R = Sp1_R_original; % first iteration 

% Interpolate demamd data on to a much finer array
MAWS_interp = 0:0.01:1400;
MV1_interp = interp1(MAWS1, MV1, MAWS_interp);
MV3_interp = interp1(MAWS3, MV3, MAWS_interp);

% convert units for easier calculations
I1 = I1 * 10^-9 * 3.154e7; % (bcf) for a year
U1original = U1original * 10^-9 * 3.154e7; % (bcf) for a year
U3original = U3original * 10^-9 * 3.154e7; % (bcf) for a year
MAWS_interp = MAWS_interp * 1e-9 * 3.154e7; % (bcf) for a year
MV1_interp = MV1_interp * 10^-6 / 43560 / 1e-9; % (milion $ / bsf) 
MV3_interp = MV3_interp * 10^-6 / 43560 / 1e-9; % (milion $ / bsf) 

% figure
% plot(MAWS1, MV1)
% 
% figure
% plot(MAWS3, MV3)
% 
% figure
% plot(MAWS_interp, MV1_interp)
% 
% figure
% plot(MAWS_interp, MV3_interp)


% Simulation loop determined watershares w/out resevoir 
for t = 1:length(Years)
    U1 = U1original;
    U3 = U3original;
    while true
        % Node 1
        u1(t) = min(U1, I1(t)); % Outflow to city
        Q1(t) = I1(t) - u1(t); % Outflow from node 1

        % Node 2
        R2(t) = 0.45 * u1(t); % Return flow from the city
        Q2(t) = Q1(t) + R2(t); % Outflow from node 2

        % Node 3
        u3(t) = min(U3, Q2(t)); % Outflow to farms
        Q3(t) = Q2(t) - u3(t); % Outflow from node 3
        
        % Environmental flow deficit test
        % if Q3(t) >= 0.5 * I1(t)
        %       0  >= 0.5 * I1(t) - Q3(t)
        Dfct(t) = 0.5 * I1(t) - Q3(t);

        if Dfct(t) < 1e-4                   
            % Calculate EW1 and EW3 for each year
            sieve1 = find(MAWS_interp <= u1(t));
            EW1(t) = (trapz(MAWS_interp(sieve1), ...
                MV1_interp(sieve1)) ...
                - MAWS_interp(sieve1(end)) ...
                * MV1_interp(sieve1(end)) ...
                - MV1_interp(sieve1(end)) ...
                * max(U1original - u1(t), 0));

            sieve3 = find(MAWS_interp <= u3(t));
            EW3(t) = (trapz(MAWS_interp(sieve3), ...
                MV3_interp(sieve3)) ...
                - MAWS_interp(sieve3(end)) ...
                * MV3_interp(sieve3(end)) ...
                - MV3_interp(sieve3(end)) ...
                * max(U3original- u3(t), 0));

            break; % Advance simulation to the next year
        else
            % Compute deficit
            Dfct(t) = 0.5 * I1(t) - Q3(t);
            % Compute reduction factors
            F1 = U1 / (U1 + U3);
            F3 = U3 / (U1 + U3);
            % Update water share targets
            U1 = max(U1 - F1 * Dfct(t), 0);
            U3 = max(U3 - F3 * Dfct(t), 0);       
        end
    end
end

% Plot Simulated flow for UIF, u1, u3, and Q3
figure 
hold on 
plot(Years, I1, 'c')
plot(Years, u1, 'r')
plot(Years, u3, 'b')
plot(Years, Q3, 'g')
hold off

% Plot EW1 and EW3
figure
hold on 
plot(Years, EW1, 'r')
plot(Years, EW3, 'b')
hold off

% Plot exceedance frequency for Simulated flow
prob_exceedance = (1:length(Years)) / (length(Years)+1);

% Plot exceedance frequency for Simulated flow for UIF, u1, u3, and Q3
rank_I1 = sort(I1, 'descend');
rank_u1 = sort(u1, 'descend');
rank_u3 = sort(u3, 'descend');
rank_Q3 = sort(Q3, 'descend');

figure
hold on 
plot(prob_exceedance(1,1:70), rank_I1, 'c')
plot(prob_exceedance(1,1:70), rank_u1, 'r')
plot(prob_exceedance(1,1:70), rank_u3, 'b')
plot(prob_exceedance(1,1:70), rank_Q3, 'g')
hold off

% Plot exceedance frequency for Simulated flow for EW1 and EW3
rank_EW1 = sort(EW1, 'descend');
rank_EW3 = sort(EW3, 'descend');

figure
hold on 
plot(prob_exceedance(1,1:70), rank_EW1, 'r')
plot(prob_exceedance(1,1:70), rank_EW3, 'b')
hold off

% % Simulation loop determined watershares w/ resevoir 
% for t = 1:length(Years)
%     U1_R = U1original;
%     U3_R = U3original; 
%     u1_R(t) = U1_R; % first iteration 
%     u3_R(t) = U3_R; % first iteration 
%     while true
%         % Node 1
%         S1_R(t+1) = S1_R(t) + I1(t) - E1_R(t) * (A1_R(t) + A1_R(t+1)) / 2 - u1_R(t) - Sp1_R(t); % Reservoir storage at end of year t        
%         R2_R(t) = 0.45 * u1_R(t);
%         Q1_R(t) = u3_R(t) - R2(t) + 0.5 * I(t); % Outflow from node 1
%         
%         % Node 2
%         R2_R(t) = 0.45 * u1_R(t); % Return flow from the city
%         Q2_R(t) = Q1_R(t) + R2_R(t) + Sp1_R(t); % Outflow from node 2
% 
%         % Node 3
%         u3_R(t) = min(U3_R, Q2_R(t)); % Outflow to farms
%         Q3_R(t) = Q2_R(t) - u3_R(t); % Outflow from node 3
% 
%         % Check for storage bounds violations:        
%         if S1_R(t+1) > S1_R_max
%             Sp1_R = S1_R(t+1) - S1_R_max;
%             if S1_R(t+1) - S1_R_max <= 1e-4 % = 0             
%                 break
%             end                 
%         end
% 
%         if S1_R(t+1) < S1_R_min
%             Dfct(t) = S1_R_min - S1(t+1);                      
%             F1 = U1_R / (U1_R + U3_R);
%             F3 = U1_R / (U1_R + U3_R);
%             u1_R = max(U1 - F1 * Dfct, 0);
%             u3_R = max(U3 - F3 * Dfct, 0);
%             if S1_R(t+1) >= -1e-4 % = S1_R_min = 0 
%                 S1_R(t+1) = S1_R_min;
%                 break
%             end  
%             S1_R(t+1) = S1_R_min;
%         end        
%     end
% end