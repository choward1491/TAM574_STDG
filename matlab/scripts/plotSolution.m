% plot solution over time of solution
% Author: C. Howard

%% specify the files to load
dir = '../../bin/';
xfile = [dir,'x_p3.txt'];
ufile = [dir,'u_p3.txt'];
qfile = [dir,'q_p3.txt'];

%% load the data
data = getSolutionData(xfile,ufile,qfile);

%% plot the final values for the concentration (q) and (u)
figure(1)
plot(data.x,data.u(end,:))
axis([data.x(1),data.x(end),-1,1])
xlabel('x','FontSize',16,'interpreter','latex')
ylabel('Pollutant Concentration, $u$','FontSize',16,'interpreter','latex')
title('$u$ vs $x$ at $t=1.0s$','FontSize',16,'interpreter','latex')

figure(2)
plot(data.x,data.q(end,:))
axis([data.x(1),data.x(end),-1,1])
xlabel('x','FontSize',16,'interpreter','latex')
ylabel('Diffusive Flux, $q$','FontSize',16,'interpreter','latex')
title('$q$ vs $x$ at $t=1.0s$','FontSize',16,'interpreter','latex')

% plot 3D visualization of space-time solution
[X,T] = meshgrid(data.x,data.t);

% concentration
figure(3)
surf(X,T,data.u)
xlabel('x','FontSize',16,'interpreter','latex')
ylabel('t','FontSize',16,'interpreter','latex')
zlabel('Pollutant Concentration, $u$','FontSize',16,'interpreter','latex')
colormap cool
shading interp
camlight

% concentration
figure(4)
surf(X,T,data.q)
xlabel('x','FontSize',16,'interpreter','latex')
ylabel('t','FontSize',16,'interpreter','latex')
zlabel('Diffusive Flux, $q$','FontSize',16,'interpreter','latex')
colormap spring
shading interp
camlight

