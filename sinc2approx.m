clc
clear all
close all

sinc2 = @(x,y)sinc(x).*sinc(y);
muAx = @(x,a,b,c)1./(1+((x-c)/a).^(2*b));
DmuAx_a = @(x,a,b,c)2*b/a^2*(x-c)*((x-c)/a)^(2*b-1)*(muAx(x,a,b,c))^2;
DmuAx_b = @(x,a,b,c)-2*(abs((x-c)/a))^(2*b)*log(abs((x-c)/a))*(muAx(x,a,b,c))^2;
DmuAx_c = @(x,a,b,c)2*b/a*((x-c)/a)^(2*b-1)*(muAx(x,a,b,c))^2;

r = 3;
step = 0.15;
[x1, x2] = meshgrid(-r:step:r, -r:step:r);
yd = sinc2(x1, x2);
trData = [x1(:) x2(:) yd(:)];            % training data set
nTrain = size(trData,1);                 % number of traning data

h1 = figure(1);
h1.Position = [8 120 560 420];
mesh(x1,x2,yd);
axis([-r r -r r -0.25 1]);
xlabel('x_{1}'); ylabel('x_{2}'); zlabel('y_{d}');
title('training data sets');

% ======== step 1. initialization of premise parameters =========
% membership function parameters for x1 fuzzification
a1 = 2*r/6*ones(1,4);           % a1 = [a11 a12 a13 a14]
b1 = 2*ones(1,4);               % b1 = [b11 b12 b13 b14]
c1 = linspace(-r+0.1,r-0.1,4);  % c1 = [c11 c12 c13 c14]

alpha1 = [a1;b1;c1];
alpha2 = alpha1;
alpha = cat(3,alpha1,alpha2);

x = -r:step:r;
h2 = figure(2);
h2.Position = [740, 400, 620, 280];
plot_MF(x, alpha, muAx);

% plotting output data with initialized parameters
A = zeros(nTrain, 48);
X = zeros(48,1);
h3 = figure(3);
h3.Position = [580 120 560 420];
Y = A*X;
g = mesh(x1,x2,reshape(Y,size(x1)));
xlabel('x_{1}'); ylabel('x_{2}'); zlabel('y');
axis([-r r -r r -0.25 1]);
pause(1);

eta = 0.2;
epoch = 2000;
Yd = trData(:,end);
E = 0.5*(Yd-Y)'*(Yd-Y);
figure(4);
ge = animatedline;
addpoints(ge,0,E);
drawnow;
axis([0 epoch/10, 0, E+1]);
xlabel('epoch'); ylabel('Error');

for e = 1:epoch
    % calculation of W
    [W, W_bar] = cal_W(alpha(:,:,1), alpha(:,:,2), trData(:,1), trData(:,2), muAx);
    
    % updating consequent parameters
    [P, Q, R] = update_param2(W, trData(:,1),trData(:,2),trData(:,3));
    Y = cal_Y(W, trData(:,1), trData(:,2), P, Q, R);
    g.ZData = reshape(Y,size(x1));
    drawnow;
    
    % Error
    E = 0.5*(Yd-Y)'*(Yd-Y);
    addpoints(ge, e, E);
    drawnow;
    % calculation of F
    F = cal_F(P,Q,R,trData(:,1),trData(:,2));
    [SW, SWF] = cal_SWF(W,F);
    
    [delta_alpha_a1, delta_alpha_a2] = update_param_a(alpha,trData(:,1),trData(:,2),eta,Yd-Y,W,F,muAx, DmuAx_a);
    [delta_alpha_b1, delta_alpha_b2] = update_param_b(alpha,trData(:,1),trData(:,2),eta,Yd-Y,W,F,muAx, DmuAx_b);
    [delta_alpha_c1, delta_alpha_c2] = update_param_c(alpha,trData(:,1),trData(:,2),eta,Yd-Y,W,F,muAx, DmuAx_c);
    alpha(1,:,1) = alpha(1,:,1) + delta_alpha_a1;
    alpha(1,:,2) = alpha(1,:,2) + delta_alpha_a2;
    alpha(2,:,1) = alpha(2,:,1) + delta_alpha_b1;
    alpha(2,:,2) = alpha(2,:,2) + delta_alpha_b2;
    alpha(3,:,1) = alpha(3,:,1) + delta_alpha_c1;
    alpha(3,:,2) = alpha(3,:,2) + delta_alpha_c2;
    
    h5 = figure(5);
    h5.Position = [740, 100, 620, 280];
    plot_MF(x, alpha, muAx);
end
