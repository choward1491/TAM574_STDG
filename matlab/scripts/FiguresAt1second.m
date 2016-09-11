%Do various polynomial figures at T = 1

close all
clear all

%Get data for p = 0
directory = '../data/';
U0 = load([directory,'u_p01.txt']);
Q0 = load([directory,'q_p01.txt']);
X0 = load([directory,'x_p01.txt']);
t = U0(:,1);
u0 = U0(:,2:end); clear U0;
q0 = Q0(:,2:end); clear Q0;

%Get data for p = 1
U1 = load([directory,'u_p11.txt']);
Q1 = load([directory,'q_p11.txt']);
X1 = load([directory,'x_p11.txt']);
u1 = U1(:,2:end); clear U1;
q1 = Q1(:,2:end); clear Q1;

%Get data for p = 2
U2 = load([directory,'u_p21.txt']);
Q2 = load([directory,'q_p21.txt']);
X2 = load([directory,'x_p21.txt']);
u2 = U2(:,2:end); clear U2;
q2 = Q2(:,2:end); clear Q2;

%Get data for p = 3
U3 = load([directory,'u_p31.txt']);
Q3 = load([directory,'q_p31.txt']);
X3 = load([directory,'x_p31.txt']);
u3 = U3(:,2:end); clear U3;
q3 = Q3(:,2:end); clear Q3;

%Get data for p = 4
U4 = load([directory,'u_p41.txt']);
Q4 = load([directory,'q_p41.txt']);
X4 = load([directory,'x_p41.txt']);
u4 = U4(:,2:end); clear U4;
q4 = Q4(:,2:end); clear Q4;

%Get data for p = 5
U5 = load([directory,'u_p51.txt']);
Q5 = load([directory,'q_p51.txt']);
X5 = load([directory,'x_p51.txt']);
u5 = U5(:,2:end); clear U5;
q5 = Q5(:,2:end); clear Q5;

%% Get the figure
for i = length(t)
    figure
    subplot(2,1,1)
    plot(X0,u0(i,:),X1,u1(i,:),'r-',X2,u2(i,:),'g.',X3,u3(i,:),'c:',X4,u4(i,:),'mx',X5,u5(i,:),'ko-')
    xlabel('Position (m)')
    ylabel('Concentration (kg/m^3)')
    title('Concentration vs Position at t = 1 second')
    legend('P = 0','P = 1','P = 2','P = 3','P = 4','P = 5',0)
    axis([0,2*pi,-.4,.4])
    subplot(2,1,2)
    plot(X0,q0(i,:),X1,q1(i,:),'r-',X2,q2(i,:),'g.',X3,q3(i,:),'c:',X4,q4(i,:),'mx',X5,q5(i,:),'ko-')
    xlabel('Position (m)')
    ylabel('Diffusive Flux (kg/(m^2s))')
    title('Diffusive Flux vs Position at t = 1 second')
    legend('P = 0','P = 1','P = 2','P = 3','P = 4','P = 5',0)
    axis([0,2*pi,-.4,.4])
    
end
