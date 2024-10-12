function plot_rr_grad_schema
%% plots schematic of the gradient ascent idea and early stopping


%% reward rate function
% scaling in x/y
rrs = [1 0.3];
rr = @(b) -sum((b(:) ./ rrs').^2);
rrgrad = @(b) -2 * (b(:)' ./ rrs.^2);
rrhess = -2 ./ rrs.^2;


%% gradient ascent properties
gradlr = 0.01;
gradini = [2 -1];
gradsteps = 15;
gradcol = [0.5 0.5 0.5];


%% plot limits and figure rotation
pxlim = [-2.5 4];
pylim = [-2 3];
prot = -10*pi/180;
protmat = [cos(prot) -sin(prot); sin(prot) cos(prot)];
xcol = [0.8 0 0];
ycol = [0 0 0.8];
graddist = 1;
gradscale = 0.1;


%% perform gradient ascent
gradp = [gradini' zeros(2, gradsteps-1)];
for n = 2:gradsteps
    gradp(:,n) = gradp(:,n-1) + gradlr * rrgrad(gradp(:,n-1))';
end


%% plot results
pbar = 1;
figure('Color', 'white');  hold on;
xlim(pxlim);  ylim(pylim);

% iso reward rates
thetas = linspace(0,2*pi,100);
isop = protmat * [(rrs(1) * cos(thetas)); (rrs(2) * sin(thetas))];
plot(isop(1,:), isop(2,:), 'k-');
plot(2 * isop(1,:), 2 * isop(2,:), 'k-');
plot(0, 0, 'k+', 'MarkerSize', 2);

% gradient ascent trajectory
gradprot = protmat * gradp;
plot(gradprot(1,:), gradprot(2,:), 'o-', 'Color', gradcol, ...
    'MarkerSize', 2, 'MarkerFaceColor', gradcol, 'MarkerEdgeColor', 'none');

% gradient and curvature in direction 1
graddir = protmat * [-(2*rrs(1)+0.2) (2*rrs(1)+graddist+0.2); 0 0];
plot(graddir(1,:), graddir(2,:), '-', 'Color', xcol);
graddir = protmat * [(gradp(1,end)-0.1) (2*rrs(1)+graddist+0.1); gradp(2,end) gradp(2,end)];
plot(graddir(1,:), graddir(2,:), '-', 'Color', xcol);
graddir = protmat * [(2*rrs(1)+graddist) (2*rrs(1)+graddist); -(2*rrs(2)+0.2) (2*rrs(2)+0.2)];
plot(graddir(1,:), graddir(2,:), '-', 'Color', xcol);
graddir = [zeros(1,50); linspace(-2*rrs(2), 2*rrs(2), 50)];
for n = 1:size(graddir,2), graddir(1,n) = rr([gradp(1,end) graddir(2,n)]); end
gradmax = max(graddir(1,:));
graddir = protmat * [(gradscale*(graddir(1,:)-gradmax)+graddist+2*rrs(1)); graddir(2,:)];
plot(graddir(1,:), graddir(2,:), 'k-');
graddir = [zeros(1,50); (gradp(2,end)+linspace(-0.3,0.3, 50))];
for n = 1:size(graddir,2), graddir(1,n) = rr([gradp(1,end) graddir(2,n)]); end
graddir = protmat * [(gradscale*(graddir(1,:)-gradmax)+graddist+2*rrs(1)); graddir(2,:)];
plot(graddir(1,:), graddir(2,:), '-', 'Color', xcol, 'LineWidth', 2);

% gradient and curvature in direction 2
graddir = protmat * [0 0; -(2*rrs(2)+0.2) (2*rrs(2)+graddist+0.2)];
plot(graddir(1,:), graddir(2,:), '-', 'Color', ycol);
graddir = protmat * [gradp(1,end) gradp(1,end); (gradp(2,end)-0.1) (2*rrs(2)+graddist+0.1)];
plot(graddir(1,:), graddir(2,:), '-', 'Color', ycol);
graddir = protmat * [-(2*rrs(1)+0.2) (2*rrs(1)+0.2); (2*rrs(2)+graddist) (2*rrs(2)+graddist)];
plot(graddir(1,:), graddir(2,:), '-', 'Color', ycol);
graddir = [linspace(-2*rrs(1), 2*rrs(1), 50); zeros(1,50)];
for n = 1:size(graddir,2), graddir(2,n) = rr([graddir(1,n) gradp(2,end)]); end
gradmax = max(graddir(2,:));
graddir = protmat * [graddir(1,:); (gradscale*(graddir(2,:)-gradmax)+graddist+2*rrs(2))];
plot(graddir(1,:), graddir(2,:), 'k-');
graddir = [(gradp(1,end)+linspace(-0.3,0.3, 50)); zeros(1,50)];
for n = 1:size(graddir,2), graddir(2,n) = rr([graddir(1,n), gradp(2,end)]); end
graddir = protmat * [graddir(1,:); (gradscale*(graddir(2,:)-gradmax)+graddist+2*rrs(2))];
plot(graddir(1,:), graddir(2,:), '-', 'Color', ycol, 'LineWidth', 2);

% set up plot
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],'XTick',[],'YTick',[]);


%% plot gradient direction plot
figure('Color', 'white');  hold on;
plot(abs(gradp(1,end)), rrhess(1), 'o', ...
    'MarkerFaceColor', xcol, 'MarkerEdgeColor', 'none', 'MarkerSize', 4);
plot(abs(gradp(2,end)), rrhess(2), 'o', ...
    'MarkerFaceColor', ycol, 'MarkerEdgeColor', 'none', 'MarkerSize', 4);
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
