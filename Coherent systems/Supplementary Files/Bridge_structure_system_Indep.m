%% Bridge structure system:

Nsamps = 100000;
KN = @(k,n) [betarnd(k, n - k + 1, Nsamps, 1), betarnd(k + 1, n - k, Nsamps, 1)];

C1 = KN(23,24); % C-box of component 1
C2 = KN(23,24); % C-box of component 2
C3 = KN(14,17); % C-box of component 3
C4 = KN(14,17); % C-box of component 4
C5 = KN(12,12); % C-box of component 5

C_L = [C1(:,1), C2(:,1), C3(:,1), C4(:,1), C5(:,1)]; % Compiling all the left side
C_R = [C1(:,2), C2(:,2), C3(:,2), C4(:,2), C5(:,2)]; % Compiling all the right side

R_L = phi_func(C_L); % Left-bounding C-box for reliability computed using structure function
R_R = phi_func(C_R); % Right-bounding C-box for reliability computed using structure function

R_L2 = Boolean_reliability(C_L); % Right-bounding C-box for reliability computed using Boolean function
R_R2 = Boolean_reliability(C_R); % Left-bounding C-box for reliability computed using Boolean function

figure; ylab = {'Probability'}; f = 25;
hold on; box on; grid on; a = 0.025; b = 0.975;
[f1,x1] = ecdf(R_L); [f2,x2] = ecdf(R_R);
stairs(x1,f1, 'r', 'linewidth', 2); stairs(x2,f2, 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'r', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'r', 'linewidth', 2, 'handlevisibility', 'off');
[f1,x1] = ecdf(R_L2); [f2,x2] = ecdf(R_R2);
stairs(x1,f1, 'b', 'linewidth', 2); stairs(x2,f2, 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([min(x1),min(x2)],[0,0], 'b', 'linewidth', 2, 'handlevisibility', 'off'); 
plot([max(x1),max(x2)],[1,1], 'b', 'linewidth', 2, 'handlevisibility', 'off');
legend('C-box obtained via Structure function', 'C-box obtained via Boolean function', 'linewidth', 2, 'location', 'northwest')
set(gca, 'Fontsize', f); xlabel('$R_{S1}$', 'Interpreter', 'latex'); ylabel(ylab); xlim([0.65,1]); title('System S1')

clc; alpha = 2.5;
confidence_interval1 = [prctile(R_L, a*100), prctile(R_R, 100-(a*100))]; disp(confidence_interval1)
confidence_interval2 = [prctile(R_L2, a*100), prctile(R_R2, 100-(a*100))]; disp(confidence_interval2)

%% Functions:

% Define the Structure function:
function [output] = phi_func(x)

output = (x(:,1) .* x(:,3)) + (x(:,2) .* x(:,4)) + (x(:,1) .* x(:,4) .* x(:,5)) + (x(:,2) .* x(:,3) .* x(:,5)) - ...
                (x(:,1) .* x(:,2) .* x(:,3) .* x(:,5)) - (x(:,1) .* x(:,3) .* x(:,4) .* x(:,5)) - ...
                (x(:,1) .* x(:,2) .* x(:,4) .* x(:,5)) - (x(:,2) .* x(:,3) .* x(:,4) .* x(:,5)) - ...
                (x(:,1) .* x(:,2) .* x(:,3) .* x(:,4)) + 2.*(x(:,1) .* x(:,2) .* x(:,3) .* x(:,4) .* x(:,5));
end

% Define the Boolean function (under independence):
function [output] = Boolean_reliability(in)
c1 = in(:,1); c2 = in(:,2); c3 = in(:,3); c4 = in(:,4); c5 = in(:,5); 

orI = @(x,y) 1 - ((1 - x).*(1 - y)); 
andI = @(x,y) x.*y;

phi1 = andI(c1,c3); phi2 = andI(c2,c4); phi3 = andI(c1,andI(c4,c5)); phi4 = andI(c2,andI(c3,c5));
output = orI((orI(phi1,phi2)), (orI(phi3,phi4)));
end