%Plot the data

clear all; close all;
flag = 0;
flag2 = 0;
GetData;
umin = min(min(u));
umax = max(max(u));
qmin = min(min(q));
qmax = max(max(q));

uexact = @(x,t) exp(-t).*sin(x);
%% The plotting
% figure
% for i = 1:length(t)
%     subplot(2,1,1)
%     plot(x,u(i,:),x,uexact(x,t(i)));
%     axis([min(x),max(x),umin,umax])
%     
%     subplot(2,1,2)
%     plot(x,q(i,:),x,-cos(x));
%    % axis([min(x),max(x),qmin,qmax])
%     pause(.01);
% end

if( flag2 == 1)
%% The Finite Difference Solution

kappa = 1;
a = 2*pi; 
tau = kappa/(a*a);
%tau = 1e-5;
x1 = linspace(0,x(end),100);
dx = x1(2)-x1(1);

t0 = linspace(0,10,100);

u0 = sin(x1)';
q0 = -kappa*cos(x1)';

IC = [u0;q0];

Dp = MatMakeExp(-3:0,length(x1),1,dx,1);
Dm = MatMakeExp(0:3,length(x1),1,dx,1);

[t,soln] = ode45(@(t,lol) TAM574(t,lol,Dp,Dm,kappa,tau,a),t,IC);
U = soln(:,1:length(x1));
Q = soln(:,length(x1)+1:end);
clear soln Dp Dm u0 q0 IC;

minu = min(min(U));
maxu = max(max(U));
minq = min(min(Q));
maxq = max(max(Q));

end
%% Plot stuff

figure
for i = 1:length(t)
        subplot(2,1,1)
        plot(x,u(i,:))
        title('Concetration vs Position')
        xlabel('Position (m)')
        ylabel('Concentration (kg/m^3)')
        axis([x(1),x(end),umin,umax])


        subplot(2,1,2)
        plot(x,q(i,:))
        title('Diffusive Flux vs Position')
        xlabel('Position (m)')
        ylabel('Diffusive Flux ( kg/ (m^2 s) )')
        axis([x(1),x(end),qmin,qmax])

    pause(.01);
end



if( flag == 1)
%% Make video
scrsz = get(0,'ScreenSize');
%fig = figure('Position',[1 scrsz(4)*.30 scrsz(3)*4/5 scrsz(4)*.7]);
fig = figure(1);
writerObj = VideoWriter('trial_p1');
open( writerObj );

    for i = 1:length(t)
        subplot(2,1,1)
        plot(x,u(i,:))
        str = ['Concetration vs Position @ t = ',num2str(t(i)),'sec' ];
        title(str)
        xlabel('Position (m)')
        ylabel('Concentration (kg/m^3)')
       % axis([x(1),x(end),umin,umax])


        subplot(2,1,2)
        plot(x,q(i,:))
        title('Diffusive Flux vs Position')
        xlabel('Position (m)')
        ylabel('Diffusive Flux ( kg/ (s*m^2) )')
        %axis([x(1),x(end),qmin,qmax])

    currFrame = getframe(fig);
    writeVideo(writerObj,currFrame);
    pause(.01);
    end
close(fig);
close(writerObj);

end