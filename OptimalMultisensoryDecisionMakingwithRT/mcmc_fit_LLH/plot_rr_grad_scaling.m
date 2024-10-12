function plot_rr_grad_scaling
%% function to show that linear re-scaling of the coordinate frame changes
%% gradient ascent

% function and scaling
f = @(x, y) - x * x - y * y;
f_grad = @(x, y) [(-2*x) (-2*y)];
%f_A_grad
A1 = [0.5 0; 0 1];
A2 = [1 0; 0 0.5];

x_ini = [-1 -1];

steps = 20;
alpha = 5e-2;


% gradient rescaling:
%
% assume that we are operating in x (vector), and we have another y = g(x) = A x
% (A is invertiable matrix). Then, it can be shown that
%
% grad_y f(g^-1(y)) = (A^-1)^T grad_x f(x),
%
% where the first is evaluated at some y, and the second at some x, such
% that y = g(x).
%
% This implies that
%
% x^n+1 = x^n + alpha grad_x f(x^n)
%
% does not map onto
%
% y^n+1 = y^n + alpha grad_y f(y^n).
%
% Instead, one finds that the latter gives
%
% x^n+1 = x^n + A^-1 (A^-1)^T grad_x f(x^n)


% perform gradient ascent simulations
X = zeros(steps, 2);
X1 = zeros(steps, 2);
X2 = zeros(steps, 2);
A1_inv = inv(A1);
A2_inv = inv(A2);
X(1,:) = x_ini;
X1(1,:) = x_ini * A1';
X2(1,:) = x_ini * A2';
for n = 2:steps
    X(n,:) = X(n-1,:) + alpha * f_grad(X(n-1,1), X(n-1,2));
    % below line re-scales the gradient such that it becomes the same as X1
    %X(n,:) = X(n-1,:) + alpha * (A1_inv * A1_inv' * f_grad(X(n-1,1), X(n-1,2))')';
    X1n1 = A1_inv * X1(n-1,:)';
    X1(n,:) = X1(n-1,:) + alpha * (A1_inv' * f_grad(X1n1(1), X1n1(2))')';
    X2n1 = A2_inv * X2(n-1,:)';
    X2(n,:) = X2(n-1,:) + alpha * (A2_inv' * f_grad(X2n1(1), X2n1(2))')';
end


% plot trajectories, X space
pbar = 1;
figure('Color', 'white'); hold on;   xlim([-1.2, 1.2]);  ylim([-1.2 1.2]); 
plot(X(:,1), X(:,2), 'ko-');
X1X = X1 * A1_inv';
plot(X1X(:, 1), X1X(:, 2), 'ro-');
X2X = X2 * A2_inv';
plot(X2X(:, 1), X2X(:, 2), 'bo-');
plot(0, 0, 'k+');
xlabel('x');
ylabel('y');
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
