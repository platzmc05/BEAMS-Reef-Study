%% Marker 32 Reef Pre_Post Restoration Ratio Analysis -  SeapHOx data 
% Michelle Platz
% 2/24/2021

close all
clc

%after running initial_look_SP_XXXX2020_M for all 5 daytime datasets
%load datasets (5 sets x 3 parameters)

    % Pre: 324
    % M_July_SDN_day = 189 points
    % M_Aug_1_SDN_day = 95 points
    % M_Aug_2_SDN_day = 40 points 
    % 
    % Post: 174
    % M_Sept_SDN_day = 122 points
    % M_Oct_SDN_day = 52 points
    
% Combine datasets
M_Pre_SDN_day = [M_July_SDN_day, M_Aug_1_SDN_day, M_Aug_2_SDN_day];
M_Pre_O2_gradient_day = [M_July_O2_gradient_day, M_Aug_1_O2_gradient_day, M_Aug_2_O2_gradient_day];
M_Pre_TA_gradient_day = [M_July_TA_gradient_day, M_Aug_1_TA_gradient_day, M_Aug_2_TA_gradient_day];

M_Post_SDN_day = [M_Sept_SDN_day, M_Oct_SDN_day];
M_Post_O2_gradient_day = [M_Sept_O2_gradient_day, M_Oct_O2_gradient_day];
M_Post_TA_gradient_day = [M_Sept_TA_gradient_day, M_Oct_TA_gradient_day];

M_Total_SDN_day = [M_Pre_SDN_day, M_Post_SDN_day];

%Plot temporal Gradient datasets 
figure
hold on; box on;

July_DOplot = plot(M_July_SDN_day, M_July_O2_gradient_day, 'b-'); %binned DO Gradient 
July_TAplot = plot(M_July_SDN_day, M_July_TA_gradient_day, 'r-'); %binned TA Gradient 

Aug_1_DOplot = plot(M_Aug_1_SDN_day, M_Aug_1_O2_gradient_day, 'b-'); %binned DO Gradient 
Aug_1_TAplot = plot(M_Aug_1_SDN_day, M_Aug_1_TA_gradient_day, 'r-'); %binned TA Gradient 

Aug_2_DOplot = plot(M_Aug_2_SDN_day, M_Aug_2_O2_gradient_day, 'b-'); %binned DO Gradient 
Aug_2_TAplot = plot(M_Aug_2_SDN_day, M_Aug_2_TA_gradient_day, 'r-'); %binned TA Gradient 

Sept_DOplot = plot(M_Sept_SDN_day, M_Sept_O2_gradient_day, 'b-'); %binned DO Gradient 
Sept_TAplot = plot(M_Sept_SDN_day, M_Sept_TA_gradient_day, 'r-'); %binned TA Gradient 

Oct_DOplot = plot(M_Oct_SDN_day, M_Oct_O2_gradient_day, 'b-'); %binned DO Gradient 
Oct_TAplot = plot(M_Oct_SDN_day, M_Oct_TA_gradient_day, 'r-'); %binned TA Gradient 

plot(M_Total_SDN_day, zeros(size(M_Total_SDN_day)), 'k'); %Zero Line

ylabel('\DeltaO_2 or \DeltaTA');
legend([July_DOplot July_TAplot], {'\DeltaO_2','\DeltaTA'}, 'location', 'northeast');
title('Marker 32 2020 Hourly Binned Concentration Gradients');


%% M Boxplots

figure
hold on; box on;
boxplot(M_Pre_TA_gradient_day)
ylabel('TA')
title('M32 Pre TA Gradient Boxplots')

figure
hold on; box on;
boxplot(M_Pre_O2_gradient_day)
ylabel('DO')
title('M32 Pre DO Gradient Boxplots')

figure
hold on; box on;
boxplot(M_Post_TA_gradient_day)
ylabel('TA')
title('M32 Post TA Gradient Boxplots')

figure
hold on; box on;
boxplot(M_Post_O2_gradient_day)
ylabel('DO')
title('M32 Post DO Gradient Boxplots')

%% Remove significant outliers  

%M_Pre:
% M_Pre_SDN_day(:,118) = [];
% M_Pre_TA_gradient_day(:,118) = [];
% M_Pre_O2_gradient_day(:,118) = [];

%M_Post: 11, 12 
% M_Post_SDN_day(:,37) = [];
% M_Post_TA_gradient_day(:,37) = [];
% M_Post_O2_gradient_day(:,37) = [];

%% Calculates NCC:NCP ratio using Geometric Mean Model II Regression 
close all 
clc

M_reg_x = [-8:1:8];
% Pre-Restoration Ratio 
[M_pre_m,M_pre_b,M_pre_r,M_pre_sm,M_pre_sb]=lsqfitgm(M_Pre_O2_gradient_day,M_Pre_TA_gradient_day)
M_Pre_Reg_Line = M_pre_m*M_reg_x + M_pre_b;
M_Pre_Ratio = M_pre_m;
M_Pre_R2 = M_pre_r;

% Post-Restoration Ratio 
[M_post_m,M_post_b,M_post_r,M_post_sm,M_post_sb]=lsqfitgm(M_Post_O2_gradient_day,M_Post_TA_gradient_day)
M_Post_Reg_Line = M_post_m*M_reg_x + M_post_b;
M_Post_Ratio = M_post_m;
M_Post_R2 = M_post_r;

M_Ratio_increase = ((M_Post_Ratio-M_Pre_Ratio)/M_Pre_Ratio)*100


%Plot Ratios
close all
%Pre-Restoration Plot
subplot(1,2,1)
hold on; 
plot(M_Pre_O2_gradient_day,M_Pre_TA_gradient_day,'bo')
plot(M_reg_x,M_Pre_Reg_Line,'r')
% xlim([-8 8])
% ylim([-5 5])
xlabel('NCP');
ylabel('NCC');
title('Marker 32 Pre-Restoration NCC:NCP Ratio');
annotation('textbox', [0.36 0.20, 0.10, 0.05], 'String', "NCC:NCP =" + M_Pre_Ratio)
annotation('textbox', [0.36, 0.15, 0.10, 0.05], 'String', "R2 =" + M_Pre_R2)

%Post-Restoration Plot
subplot(1,2,2)
hold on; 
plot(M_Post_O2_gradient_day,M_Post_TA_gradient_day,'bo')
plot(M_reg_x,M_Post_Reg_Line,'r')
% xlim([-8 8])
% ylim([-5 5])
xlabel('NCP');
ylabel('NCC');
title('Marker 32 Post-Restoration NCC:NCP Ratio');
annotation('textbox', [0.8 0.20, 0.10, 0.05], 'String', "NCC:NCP =" + M_Post_Ratio)
annotation('textbox', [0.8, 0.15, 0.10, 0.05], 'String', "R2 =" + M_Post_R2)
annotation('textbox', [0.8, 0.25, 0.10, 0.05], 'String', "Percent Ratio Increase =" + M_Ratio_increase +'%')

%% ***** Cudjoe Reef Pre_Post Restoration Ratio Analysis -  SeapHOx data 
% Michelle Platz
% 2/25/2021
% updated 3/16/2021 

close all
clc

%after running initial_look_SP_XXXX2020_U for all 5 daytime datasets
%load datasets (5 sets x 3 parameters)

%split Aug_2 into pre- and post- restoration -- 1-66 = pre, 67-125 = post

U_Aug_2_SDN_day_pre = U_Aug_2_SDN_day(1:66);
U_Aug_2_O2_gradient_day_pre = U_Aug_2_O2_gradient_day(1:66);
U_Aug_2_TA_gradient_day_pre = U_Aug_2_TA_gradient_day(1:66);

U_Aug_2_SDN_day_post = U_Aug_2_SDN_day(67:125);
U_Aug_2_O2_gradient_day_post = U_Aug_2_O2_gradient_day(67:125);
U_Aug_2_TA_gradient_day_post = U_Aug_2_TA_gradient_day(67:125);

    % Pre: 202
    % U_July_SDN_day = 37 points
    % U_Aug_1_SDN_day = 99 points
    % U_Aug_2_SDN_day_pre = 66 points 
    % 
    % Post: 290
    % U_Aug_2_SDN_day_post = 59 points
    % U_Sept_SDN_day = 129 points
    % U_Oct_SDN_day = 102 points
    
% Combine datasets
% leave out U_July 
U_Pre_SDN_day = [U_July_SDN_day, U_Aug_1_SDN_day, U_Aug_2_SDN_day_pre];
U_Pre_O2_gradient_day = [U_July_O2_gradient_day, U_Aug_1_O2_gradient_day, U_Aug_2_O2_gradient_day_pre];
U_Pre_TA_gradient_day = [U_July_TA_gradient_day, U_Aug_1_TA_gradient_day, U_Aug_2_TA_gradient_day_pre];

U_Post_SDN_day = [U_Aug_2_SDN_day_post, U_Sept_SDN_day, U_Oct_SDN_day1, U_Oct_SDN_day2];
U_Post_O2_gradient_day = [U_Aug_2_O2_gradient_day_post, U_Sept_O2_gradient_day, U_Oct_O2_gradient_day1, U_Oct_O2_gradient_day2];
U_Post_TA_gradient_day = [U_Aug_2_TA_gradient_day_post, U_Sept_TA_gradient_day, U_Oct_TA_gradient_day1, U_Oct_TA_gradient_day2];

U_Total_SDN_day = [U_Pre_SDN_day, U_Post_SDN_day];

%Plot temporal Gradient datasets - day only 
figure
hold on; box on;

July_DOplot = plot(U_July_SDN_day, U_July_O2_gradient_day, 'b-'); %binned DO Gradient 
July_TAplot = plot(U_July_SDN_day, U_July_TA_gradient_day, 'r-'); %binned TA Gradient 

Aug_1_DOplot = plot(U_Aug_1_SDN_day, U_Aug_1_O2_gradient_day, 'b-'); %binned DO Gradient 
Aug_1_TAplot = plot(U_Aug_1_SDN_day, U_Aug_1_TA_gradient_day, 'r-'); %binned TA Gradient 

Aug_2_DOplot = plot(U_Aug_2_SDN_day, U_Aug_2_O2_gradient_day, 'b-'); %binned DO Gradient 
Aug_2_TAplot = plot(U_Aug_2_SDN_day, U_Aug_2_TA_gradient_day, 'r-'); %binned TA Gradient 

Sept_DOplot = plot(U_Sept_SDN_day, U_Sept_O2_gradient_day, 'b-'); %binned DO Gradient 
Sept_TAplot = plot(U_Sept_SDN_day, U_Sept_TA_gradient_day, 'r-'); %binned TA Gradient 

Oct_DOplot1 = plot(U_Oct_SDN_day1, U_Oct_O2_gradient_day1, 'b-'); %binned DO Gradient 
Oct_TAplot1 = plot(U_Oct_SDN_day1, U_Oct_TA_gradient_day1, 'r-'); %binned TA Gradient 

Oct_DOplot2 = plot(U_Oct_SDN_day2, U_Oct_O2_gradient_day2, 'b-'); %binned DO Gradient 
Oct_TAplot2 = plot(U_Oct_SDN_day2, U_Oct_TA_gradient_day2, 'r-'); %binned TA Gradient 

plot(U_Total_SDN_day, zeros(size(U_Total_SDN_day)), 'k'); %Zero Line

ylabel('\DeltaO_2 or \DeltaTA');
legend([July_DOplot July_TAplot], {'\DeltaO_2','\DeltaTA'}, 'location', 'northeast');
title('Cudjoe 2020 Hourly Binned Concentration Gradients');


%% U Boxplots
close all 

figure
hold on; box on;
boxplot(U_Pre_TA_gradient_day)
ylabel('TA')
title('U Pre TA Gradient Boxplots')

figure
hold on; box on;
boxplot(U_Pre_O2_gradient_day)
ylabel('DO')
title('U Pre DO Gradient Boxplots')

figure
hold on; box on;
boxplot(U_Post_TA_gradient_day)
ylabel('TA')
title('U Post TA Gradient Boxplots')

figure
hold on; box on;
boxplot(U_Post_O2_gradient_day)
ylabel('DO')
title('U Post DO Gradient Boxplots')

%% Remove significant outliers  
% 
% %U_Pre: 7
% U_Pre_SDN_day(:,7) = [];
% U_Pre_TA_gradient_day(:,7) = [];
% U_Pre_O2_gradient_day(:,7) = [];
% % 
% % %U_Post: 11, 12 
% U_Post_SDN_day(:,11) = [];
% U_Post_TA_gradient_day(:,11) = [];
% U_Post_O2_gradient_day(:,11) = [];


%% Calculates NCC:NCP ratio using Geometric Mean Model II Regression 
close all 
clc

U_reg_x = [-8:1:8];
% Pre-Restoration Ratio 
[U_pre_m,U_pre_b,U_pre_r,U_pre_sm,U_pre_sb]=lsqfitgm(U_Pre_O2_gradient_day,U_Pre_TA_gradient_day)
U_Pre_Reg_Line = U_pre_m*U_reg_x + U_pre_b;
U_Pre_Ratio = U_pre_m;
U_Pre_R2 = U_pre_r;

% Post-Restoration Ratio 
[U_post_m,U_post_b,U_post_r,U_post_sm,U_post_sb]=lsqfitgm(U_Post_O2_gradient_day,U_Post_TA_gradient_day)
U_Post_Reg_Line = U_post_m*U_reg_x + U_post_b;
U_Post_Ratio = U_post_m;
U_Post_R2 = U_post_r;

%Calculate % increase in ratio
U_Ratio_increase = ((U_Post_Ratio-U_Pre_Ratio)/U_Pre_Ratio)*100

% Plot Ratios
clc
%Pre-Restoration Plot
subplot(1,2,1)
hold on; 
plot(U_Pre_O2_gradient_day,U_Pre_TA_gradient_day,'bo')
plot(U_reg_x,U_Pre_Reg_Line,'r')
% xlim([-5 5])
% ylim([-2 2])
xlabel('NCP');
ylabel('NCC');
title('Cudjoe Pre-Restoration NCC:NCP Ratio');
annotation('textbox', [0.36 0.20, 0.10, 0.05], 'String', "NCC:NCP =" + U_Pre_Ratio)
annotation('textbox', [0.36, 0.15, 0.10, 0.05], 'String', "R2 =" + U_Pre_R2)

%Post-Restoration Plot
subplot(1,2,2)
hold on; 
plot(U_Post_O2_gradient_day,U_Post_TA_gradient_day,'bo')
plot(U_reg_x,U_Post_Reg_Line,'r')
% xlim([-5 5])
% ylim([-2 2])
xlabel('NCP');
ylabel('NCC');
title('Cudjoe Post-Restoration NCC:NCP Ratio');
annotation('textbox', [0.8 0.20, 0.10, 0.05], 'String', "NCC:NCP =" + U_Post_Ratio)
annotation('textbox', [0.8, 0.15, 0.10, 0.05], 'String', "R2 =" + U_Post_R2)
annotation('textbox', [0.8, 0.25, 0.10, 0.05], 'String', "Percent Ratio Increase =" + U_Ratio_increase +'%')







