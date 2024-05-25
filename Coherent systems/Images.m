%% Define the k-out-of-n C-box:

Nsamps = 100000;
KN = @(k,n) [betarnd(k, n - k + 1, Nsamps, 1), betarnd(k + 1, n - k, Nsamps, 1)];

%% Illustration

samp = KN(6,10);

figure; ylab = {'Probability'}; f =18; %
hold on; box on; grid on; a = 0.025; b = 0.975;
perctile = prctile(samp, [a*100, 100-(a*100)]); conf_int = [perctile(1,1), perctile(2,2)];
[f1,x1] = ecdf(samp(:,1)); [f2,x2] = ecdf(samp(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2); 
plot([0, conf_int(1)], [a, a], 'k--', 'linewidth', 2); plot([0, conf_int(2)], [b, b], 'm--', 'linewidth', 2);
plot([conf_int(1), conf_int(1)], [0, a], 'k--', 'linewidth', 2, 'handlevisibility', 'off'); plot([conf_int(2), conf_int(2)], [0, b], 'm--', 'linewidth', 2, 'handlevisibility', 'off');
plot([conf_int(1), conf_int(2)], [0.004, 0.004], 'g', 'linewidth', 4, 'handlevisibility', 'off');
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('$\theta$', 'Interpreter', 'latex'); ylabel(ylab,'FontName', 'Times'); xlim([0,1]);
legend('beta(6, 5)', 'beta(7, 4)', '\alpha = 0.025', '\beta = 0.975', 'linewidth', 2, 'location', 'southeast','FontName', 'Times')

%% Example: Omega component reliability:

C = KN(2, 60);

figure; ylab = {'Probability'}; f =18;
hold on; box on; grid on; a = 0.05; b = 0.95;
perctile = prctile(C, [a*100, 100-(a*100)]); conf_int = [perctile(1,1), perctile(2,2)];
[f1,x1] = ecdf(C(:,1)); [f2,x2] = ecdf(C(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2); 
plot([0, conf_int(1)], [a, a], 'k--', 'linewidth', 2); plot([0, conf_int(2)], [b, b], 'm--', 'linewidth', 2);
plot([conf_int(1), conf_int(1)], [0, a], 'k--', 'linewidth', 2, 'handlevisibility', 'off'); plot([conf_int(2), conf_int(2)], [0, b], 'm--', 'linewidth', 2, 'handlevisibility', 'off');
plot([conf_int(1), conf_int(2)], [0.004, 0.004], 'g', 'linewidth', 4, 'handlevisibility', 'off');
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('$p_{f}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.2]); title('Product \Omega', 'FontName', 'Times')
legend('beta(2, 59)', 'beta(3, 58)', '\alpha = 0.05', '\beta = 0.95', 'linewidth', 2, 'location', 'southeast', 'FontName', 'Times')

%% Example: Series / Parallel system:

C1 = KN(23, 24); % Component #1
C2 = KN(14, 17); % Component #2

ser(:,:,1) = csvread("table_R_indep_series.csv", 1, 1); ser(:,:,2) = csvread("table_R_frechet_series.csv", 1, 1); 
par(:,:,1) = csvread("table_R_indep_par.csv", 1, 1); par(:,:,2) = csvread("table_R_frechet_par.csv", 1, 1); 

figure; ylab = {'Probability'}; f =18;
subplot(2,2,1)
hold on; box on; grid on;
[f1,x1] = ecdf(C1(:,1)); [f2,x2] = ecdf(C1(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('$p_{1}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0.6,1]); title('Component 1', 'FontName', 'Times')

subplot(2,2,2)
hold on; box on; grid on;
[f1,x1] = ecdf(C2(:,1)); [f2,x2] = ecdf(C2(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('$p_{2}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0.3,1]); title('Component 2', 'FontName', 'Times')

subplot(2,2,3)
hold on; box on; grid on;
[f1,x1] = ecdf(ser(:,1,1)); [f2,x2] = ecdf(ser(:,2,1));
stairs(x1,f1, 'g', 'linewidth', 2); stairs(x2,f2, 'g', 'linewidth', 2, 'handlevisibility', 'off');
plot([min(x1),min(x2)],[0,0], 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'g', 'linewidth', 2, 'handlevisibility', 'off');
[f1,x1] = ecdf(ser(:,1,2)); [f2,x2] = ecdf(ser(:,2,2));
stairs(x1,f1, 'b', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'b', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('$R_{series}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,1]); ylim([0,1]); title('Series system', 'FontName', 'Times')

subplot(2,2,4)
hold on; box on; grid on;
[f1,x1] = ecdf(par(:,1,1)); [f2,x2] = ecdf(par(:,2,1));
stairs(x1,f1, 'g', 'linewidth', 2); stairs(x2,f2, 'g', 'linewidth', 2, 'handlevisibility', 'off');
plot([min(x1),min(x2)],[0,0], 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'g', 'linewidth', 2, 'handlevisibility', 'off');
[f1,x1] = ecdf(par(:,1,2)); [f2,x2] = ecdf(par(:,2,2));
stairs(x1,f1, 'b', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'b', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('$R_{parallel}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,1]); ylim([0,1]); title('Parallel system', 'FontName', 'Times')
legend('Independence', 'Uncertain dependencies', 'linewidth', 2, 'FontName', 'Times')

%% Case Study 1: Pressurised Tank system:

C1 = KN(0, 2000);  % Pressure tank T
C2 = KN(3, 500);   % Relay K2
C3 = KN(1, 150);   % Pressure switch S
C4 = KN(0, 460);   % On-off switch S1
C5 = KN(7, 10000); % Relay K1
C6 = KN(0, 380);   % Timer relay R

figure; ylab = {'Probability'}; f = 18;
subplot(2,3,1)
hold on; box on; grid on;
[f1,x1] = ecdf(C1(:,1)); [f2,x2] = ecdf(C1(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{1}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.0034]); title('Pressure tank T', 'FontName', 'Times')

subplot(2,3,2)
hold on; box on; grid on;
[f1,x1] = ecdf(C2(:,1)); [f2,x2] = ecdf(C2(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{2}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.02]); title('Relay K2', 'FontName', 'Times')

subplot(2,3,3)
hold on; box on; grid on;
[f1,x1] = ecdf(C3(:,1)); [f2,x2] = ecdf(C3(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{3}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.05]); title('Pressure switch S', 'FontName', 'Times')

subplot(2,3,4)
hold on; box on; grid on;
[f1,x1] = ecdf(C4(:,1)); [f2,x2] = ecdf(C4(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{4}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.013]); title('On-off switch S1', 'FontName', 'Times')

subplot(2,3,5)
hold on; box on; grid on;
[f1,x1] = ecdf(C5(:,1)); [f2,x2] = ecdf(C5(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{5}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.0018]); title('Relay K1', 'FontName', 'Times')

subplot(2,3,6)
hold on; box on; grid on;
[f1,x1] = ecdf(C6(:,1)); [f2,x2] = ecdf(C6(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{6}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.015]); title('Timer relay R', 'FontName', 'Times')

E1(:,:,1) = csvread("E1_cbox_indep.csv", 1, 1); E1(:,:,2) = csvread("E1_cbox_frechet.csv", 1, 1); 
E2(:,:,1) = csvread("E2_cbox_indep.csv", 1, 1); E2(:,:,2) = csvread("E2_cbox_frechet.csv", 1, 1); 
E3(:,:,1) = csvread("E3_cbox_indep.csv", 1, 1); E3(:,:,2) = csvread("E3_cbox_frechet.csv", 1, 1); 
E4(:,:,1) = csvread("E4_cbox_indep.csv", 1, 1); E4(:,:,2) = csvread("E4_cbox_frechet.csv", 1, 1); 
E5(:,:,1) = csvread("E5_cbox_indep.csv", 1, 1); E5(:,:,2) = csvread("E5_cbox_frechet.csv", 1, 1); 

figure; ylab = {'Probability'}; f = 18;
subplot(2,3,1)
hold on; box on; grid on;
[f1,x1] = ecdf(E1(:,1,1)); [f2,x2] = ecdf(E1(:,2,1));
stairs(x1,f1, 'g', 'linewidth', 2); stairs(x2,f2, 'g', 'linewidth', 2, 'handlevisibility', 'off');
plot([min(x1),min(x2)],[0,0], 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'g', 'linewidth', 2, 'handlevisibility', 'off');
[f1,x1] = ecdf(E1(:,1,2)); [f2,x2] = ecdf(E1(:,2,2));
stairs(x1,f1, 'b', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'b', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{E1}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.04]); ylim([0,1]); title('E1', 'FontName', 'Times')

subplot(2,3,2)
hold on; box on; grid on;
[f1,x1] = ecdf(E2(:,1,1)); [f2,x2] = ecdf(E2(:,2,1));
stairs(x1,f1, 'g', 'linewidth', 2); stairs(x2,f2, 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'g', 'linewidth', 2, 'handlevisibility', 'off');
[f1,x1] = ecdf(E2(:,1,2)); [f2,x2] = ecdf(E2(:,2,2));
stairs(x1,f1, 'b', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'b', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{E2}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.04]); ylim([0,1]); title('E2', 'FontName', 'Times')

subplot(2,3,3)
hold on; box on; grid on;
[f1,x1] = ecdf(E3(:,1,1)); [f2,x2] = ecdf(E3(:,2,1));
stairs(x1,f1, 'g', 'linewidth', 2); stairs(x2,f2, 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'g', 'linewidth', 2, 'handlevisibility', 'off');
[f1,x1] = ecdf(E3(:,1,2)); [f2,x2] = ecdf(E3(:,2,2));
stairs(x1,f1, 'b', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'b', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{E3}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.02]); ylim([0,1]); title('E3', 'FontName', 'Times')

subplot(2,3,4)
hold on; box on; grid on;
[f1,x1] = ecdf(E4(:,1,1)); [f2,x2] = ecdf(E4(:,2,1));
stairs(x1,f1, 'g', 'linewidth', 2); stairs(x2,f2, 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'g', 'linewidth', 2, 'handlevisibility', 'off');
[f1,x1] = ecdf(E4(:,1,2)); [f2,x2] = ecdf(E4(:,2,2));
stairs(x1,f1, 'b', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'b', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{E4}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.03]); ylim([0,1]); title('E4', 'FontName', 'Times')

subplot(2,3,5)
hold on; box on; grid on;
[f1,x1] = ecdf(E5(:,1,1)); [f2,x2] = ecdf(E5(:,2,1));
stairs(x1,f1, 'g', 'linewidth', 2); stairs(x2,f2, 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'g', 'linewidth', 2, 'handlevisibility', 'off');
[f1,x1] = ecdf(E5(:,1,2)); [f2,x2] = ecdf(E5(:,2,2));
stairs(x1,f1, 'b', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'b', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{E5}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.03]); ylim([0,1]); title('E5', 'FontName', 'Times')
legend('Independence', 'Uncertain dependencies', 'linewidth', 2, 'FontName', 'Times')

%% Case Study 2: TRIGA Nuclear Research Reactor:

C1 = KN(1, 42); % Inlet pipeline
C2 = KN(8, 42); % Valve V-3
C3 = 3.49E-05;  % Operator of Valve V-3
C4 = KN(8, 42); % Valve V-4
C5 = 3.49E-05;  % Operator of Valve V-4

figure; ylab = {'Probability'}; f = 18;
subplot(2,3,1)
hold on; box on; grid on;
[f1,x1] = ecdf(C1(:,1)); [f2,x2] = ecdf(C1(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{1}$ $[yr^{-1}]$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.35]); title('Inlet pipeline', 'FontName', 'Times')

subplot(2,3,2)
hold on; box on; grid on;
[f1,x1] = ecdf(C2(:,1)); [f2,x2] = ecdf(C2(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([0,min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{2}$ $[yr^{-1}]$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.6]); title('Valve V-3', 'FontName', 'Times')

subplot(2,3,3)
hold on; box on; grid on;
[f1,x1] = ecdf(C3);
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([0,min(x1)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1), 1],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{3}$ $[yr^{-1}]$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([2E-05,5E-05]); title('Operator of Valve V-3', 'FontName', 'Times')

subplot(2,3,4)
hold on; box on; grid on;
[f1,x1] = ecdf(C4(:,1)); [f2,x2] = ecdf(C4(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{4}$ $[yr^{-1}]$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.6]); title('Valve V-4', 'FontName', 'Times')

subplot(2,3,5)
hold on; box on; grid on;
[f1,x1] = ecdf(C5);
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([0,min(x1)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1), 1],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{5}$ $[yr^{-1}]$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([2E-05,5E-05]); title('Operator of Valve V4', 'FontName', 'Times')

E1(:,:,1) = csvread("E1_cbox_indep_TRIGA.csv", 1, 1); E1(:,:,2) = csvread("E1_cbox_frechet_TRIGA.csv", 1, 1); 
E2(:,:,1) = csvread("E2_cbox_indep_TRIGA.csv", 1, 1); E2(:,:,2) = csvread("E2_cbox_frechet_TRIGA.csv", 1, 1); 
E3(:,:,1) = csvread("E3_cbox_indep_TRIGA.csv", 1, 1); E3(:,:,2) = csvread("E3_cbox_frechet_TRIGA.csv", 1, 1); 
E4(:,:,1) = csvread("E4_cbox_indep_TRIGA.csv", 1, 1); E4(:,:,2) = csvread("E4_cbox_frechet_TRIGA.csv", 1, 1); 

figure; ylab = {'Probability'}; f = 18;
subplot(2,2,1)
hold on; box on; grid on;
[f1,x1] = ecdf(E1(:,1,1)); [f2,x2] = ecdf(E1(:,2,1));
stairs(x1,f1, 'g', 'linewidth', 2); stairs(x2,f2, 'g', 'linewidth', 2, 'handlevisibility', 'off');
plot([min(x1),min(x2)],[0,0], 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'g', 'linewidth', 2, 'handlevisibility', 'off');
[f1,x1] = ecdf(E1(:,1,2)); [f2,x2] = ecdf(E1(:,2,2));
stairs(x1,f1, 'b', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'b', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{E1}$ $[yr^{-1}]$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.2]); ylim([0,1]); title('E1', 'FontName', 'Times')

subplot(2,2,2)
hold on; box on; grid on;
[f1,x1] = ecdf(E2(:,1,1)); [f2,x2] = ecdf(E2(:,2,1));
stairs(x1,f1, 'g', 'linewidth', 2); stairs(x2,f2, 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'g', 'linewidth', 2, 'handlevisibility', 'off');
[f1,x1] = ecdf(E2(:,1,2)); [f2,x2] = ecdf(E2(:,2,2));
stairs(x1,f1, 'b', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'b', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{E2}$ $[yr^{-1}]$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,1]); ylim([0,1]); title('E2', 'FontName', 'Times')

subplot(2,2,3)
hold on; box on; grid on;
[f1,x1] = ecdf(E3(:,1,1)); [f2,x2] = ecdf(E3(:,2,1));
stairs(x1,f1, 'g', 'linewidth', 2); stairs(x2,f2, 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'g', 'linewidth', 2, 'handlevisibility', 'off');
[f1,x1] = ecdf(E3(:,1,2)); [f2,x2] = ecdf(E3(:,2,2));
stairs(x1,f1, 'b', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'b', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{E3}$ $[yr^{-1}]$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.4]); ylim([0,1]); title('E3', 'FontName', 'Times')

subplot(2,2,4)
hold on; box on; grid on;
[f1,x1] = ecdf(E4(:,1,1)); [f2,x2] = ecdf(E4(:,2,1));
stairs(x1,f1, 'g', 'linewidth', 2); stairs(x2,f2, 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'g', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'g', 'linewidth', 2, 'handlevisibility', 'off');
[f1,x1] = ecdf(E4(:,1,2)); [f2,x2] = ecdf(E4(:,2,2));
stairs(x1,f1, 'b', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'b', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('${p_{f}}_{E4}$ $[yr^{-1}]$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0,0.4]); ylim([0,1]); title('E4', 'FontName', 'Times')
legend('Independence', 'Uncertain dependencies', 'linewidth', 2, 'FontName', 'Times')

%% Case Study 3: Bridge system:

C1 = KN(23, 24); % Component #1
C2 = KN(23, 24); % Component #2
C3 = KN(14, 17); % Component #3
C4 = KN(14, 17); % Component #4
C5 = KN(12, 12); % Component #5

figure; ylab = {'Probability'}; f = 18;
subplot(2,3,1)
hold on; box on; grid on;
[f1,x1] = ecdf(C1(:,1)); [f2,x2] = ecdf(C1(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('$p_{1}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0.6,1]); title('Component 1', 'FontName', 'Times')

subplot(2,3,2)
hold on; box on; grid on;
[f1,x1] = ecdf(C2(:,1)); [f2,x2] = ecdf(C2(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('$p_{2}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0.6,1]); title('Component 2', 'FontName', 'Times')

subplot(2,3,3)
hold on; box on; grid on;
[f1,x1] = ecdf(C3(:,1)); [f2,x2] = ecdf(C3(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('$p_{3}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0.3,1]); title('Component 3', 'FontName', 'Times')

subplot(2,3,4)
hold on; box on; grid on;
[f1,x1] = ecdf(C4(:,1)); [f2,x2] = ecdf(C4(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('$p_{4}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0.35,1]); title('Component 4', 'FontName', 'Times')

subplot(2,3,5)
hold on; box on; grid on;
[f1,x1] = ecdf(C5(:,1)); [f2,x2] = ecdf(C5(:,2));
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
set(gca, 'Fontsize', f); xlabel('$p_{5}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0.4,1]); title('Component 5', 'FontName', 'Times')

% Monte Carlo approximation of Bridge system structure reliability under independence:
C_L = [C1(:,1), C2(:,1), C3(:,1), C4(:,1), C5(:,1)]; % Compiling all the left side
C_R = [C1(:,2), C2(:,2), C3(:,2), C4(:,2), C5(:,2)]; % Compiling all the right side

% Define the structure function:
phi_func = @(x) (x(:,1) .* x(:,3)) + (x(:,2) .* x(:,4)) + (x(:,1) .* x(:,4) .* x(:,5)) + (x(:,2) .* x(:,3) .* x(:,5)) - ...
                (x(:,1) .* x(:,2) .* x(:,3) .* x(:,5)) - (x(:,1) .* x(:,3) .* x(:,4) .* x(:,5)) - ...
                (x(:,1) .* x(:,2) .* x(:,4) .* x(:,5)) - (x(:,2) .* x(:,3) .* x(:,4) .* x(:,5)) - ...
                (x(:,1) .* x(:,2) .* x(:,3) .* x(:,4)) + 2.*(x(:,1) .* x(:,2) .* x(:,3) .* x(:,4) .* x(:,5));

R_L = phi_func(C_L); R_R = phi_func(C_R);
R = csvread("R_case1.csv", 1, 1);

figure;
hold on; box on; grid on; f = 18; col = {'b', 'm'}; ylab = {'Probability'};
[f1,x1] = ecdf(R_L); [f2,x2] = ecdf(R_R);
stairs(x1, f1, 'color', col{1}, 'linewidth', 2); stairs(x2, f2, 'color', col{1}, 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)], [0,0], 'color', col{1}, 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)], [1,1], 'color', col{1}, 'linewidth', 2, 'handlevisibility', 'off');

[f1,x1] = ecdf(R(:,1)); [f2,x2] = ecdf(R(:,2));
stairs(x1, f1, 'color', col{2}, 'linewidth', 2); stairs(x2, f2, 'color', col{2}, 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)], [0,0], 'color', col{2}, 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)], [1,1], 'color', col{2}, 'linewidth', 2, 'handlevisibility', 'off');

set(gca, 'Fontsize', f); xlabel('$R_{S1}$', 'Interpreter', 'latex'); ylabel(ylab, 'FontName', 'Times'); xlim([0.5,1]); title('System S1', 'FontName', 'Times')
legend('Structure function', 'Boolean expression', 'linewidth', 2, 'location', 'northwest', 'FontName', 'Times')

% Compute confidence interval for Monte Carlo method
a = 0.025; b = 0.975;
confidence_interval_MC = [prctile(R_L, a*100), prctile(R_R, 100-(a*100))];
confidence_intervalBoolean = [prctile(R(:,1), a*100), prctile(R(:,2), 100-(a*100))];

%%
save('Bridge_structure_system')