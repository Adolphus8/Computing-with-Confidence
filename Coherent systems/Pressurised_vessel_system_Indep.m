%% Pressure vessel system: 
% Comparing structure function solution to that of Boolean expression
%%

Nsamps = 100000;
KN = @(k,n) [betarnd(k, n - k + 1, Nsamps, 1), betarnd(k + 1, n - k, Nsamps, 1)];

C1 = KN(2, 2000); % Pressure tank T
C2 = KN(3, 500); % Relay K2
C3 = KN(1, 150); % Pressure switch S
C4 = KN(0, 460); % On-off switch S1
C5 = KN(7, 10000); % Relay K1
C6 = KN(0, 380); % Timer relay R

C_L = [C1(:,1), C2(:,1), C3(:,1), C4(:,1), C5(:,1), C6(:,1)]; % Compiling all the left side
C_R = [C1(:,2), C2(:,2), C3(:,2), C4(:,2), C5(:,2), C6(:,2)]; % Compiling all the right side

R_L = phi_func(C_L); % Left-bounding C-box for reliability computed using structure function
R_R = phi_func(C_R); % Right-bounding C-box for reliability computed using structure function

R_L2 = Boolean_reliability(C_L); % Right-bounding C-box for reliability computed using Boolean function
R_R2 = Boolean_reliability(C_R); % Left-bounding C-box for reliability computed using Boolean function

figure; ylab = {'Probability'}; f = 25;
hold on; box on; grid on; a = 0.025; b = 0.975;
[f1,x1] = ecdf(R_L); [f2,x2] = ecdf(R_R);
stairs(x1,f1, 'r', 'linewidth', 3); stairs(x2,f2, 'r', 'linewidth', 3, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 3, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 3, 'handlevisibility', 'off');
[f1,x1] = ecdf(R_L2); [f2,x2] = ecdf(R_R2);
stairs(x1,f1, 'b', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'b', 'linewidth', 2, 'handlevisibility', 'off');
legend('C-box obtained via Structure function', 'C-box obtained via Boolean function', 'linewidth', 2, 'location', 'southeast')
set(gca, 'Fontsize', f); xlabel('$p_{f}$', 'Interpreter', 'latex'); ylabel(ylab); xlim([0, 0.03]); title('Pressure Vessel System')

clc; alpha = 2.5;
confidence_interval1 = [prctile(R_L, a*100), prctile(R_R, 100-(a*100))]; disp(confidence_interval1)
confidence_interval2 = [prctile(R_L2, a*100), prctile(R_R2, 100-(a*100))]; disp(confidence_interval2)

%% Functions:

% Define the Structure function:
function [output] = phi_func(x)

output = x(:,1) + x(:,2) + (x(:,3).*x(:,4)) + (x(:,3).*x(:,5)) + (x(:,3).*x(:,6)) - (x(:,3).*x(:,5).*x(:,6)) - ...
         (x(:,3).*x(:,4).*x(:,6)) - (x(:,3).*x(:,4).*x(:,5)) + (x(:,3).*x(:,4).*x(:,5).*x(:,6)) - (x(:,2).*x(:,3).*x(:,4)) - ...
         (x(:,2).*x(:,3).*x(:,5)) + (x(:,2).*x(:,3).*x(:,5).*x(:,6)) + (x(:,2).*x(:,3).*x(:,4).*x(:,6)) + ...
         (x(:,2).*x(:,3).*x(:,4).*x(:,5)) - (x(:,2).*x(:,3).*x(:,4).*x(:,5).*x(:,6)) - (x(:,1).*x(:,2)) - (x(:,1).*x(:,3).*x(:,4)) - ...
         (x(:,1).*x(:,3).*x(:,5)) - (x(:,1).*x(:,3).*x(:,6)) + (x(:,1).*x(:,3).*x(:,5).*x(:,6)) + (x(:,1).*x(:,3).*x(:,4).*x(:,6)) + ...
         (x(:,1).*x(:,3).*x(:,4).*x(:,5)) - (x(:,1).*x(:,2).*x(:,3).*x(:,5).*x(:,6)) - (x(:,1).*x(:,2).*x(:,3).*x(:,4).*x(:,6)) - ...
         (x(:,1).*x(:,2).*x(:,3).*x(:,4).*x(:,5)) + (x(:,1).*x(:,2).*x(:,3).*x(:,4).*x(:,5).*x(:,6));
end

% Define the Boolean function (under independence):
function [output] = Boolean_reliability(in)
c1 = in(:,1); c2 = in(:,2); c3 = in(:,3); c4 = in(:,4); c5 = in(:,5); c6 = in(:,6); 

orI = @(x,y) 1 - ((1 - x).*(1 - y)); 
andI = @(x,y) x.*y;

E5 = orI(c6, c5); 
E4 = orI(E5, c4);
E3 = andI(E4, c3);
E2 = orI(E3, c2);
output = orI(E2,c1);
end