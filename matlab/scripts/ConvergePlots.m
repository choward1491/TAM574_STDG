%Do various polynomial figures at T = 1

close all
clear all

%Get data for p = 0
directory = '../data/';
U01 = load([directory,'u_p0_c1.txt']);
U02 = load([directory,'u_p0_c2.txt']);
U04 = load([directory,'u_p0_c4.txt']);
U08 = load([directory,'u_p0_c8.txt']);
Q01 = load([directory,'q_p0_c1.txt']);
Q02 = load([directory,'q_p0_c2.txt']);
Q04 = load([directory,'q_p0_c4.txt']);
Q08 = load([directory,'q_p0_c8.txt']);
X01 = load([directory,'x_p0_c1.txt']);
X02 = load([directory,'x_p0_c2.txt']);
X04 = load([directory,'x_p0_c4.txt']);
X08 = load([directory,'x_p0_c8.txt']);


%Get data for p = 1
U11 = load([directory,'u_p1_c1.txt']);
U12 = load([directory,'u_p1_c2.txt']);
U14 = load([directory,'u_p1_c4.txt']);
U18 = load([directory,'u_p1_c8.txt']);
Q11 = load([directory,'q_p1_c1.txt']);
Q12 = load([directory,'q_p1_c2.txt']);
Q14 = load([directory,'q_p1_c4.txt']);
Q18 = load([directory,'q_p1_c8.txt']);
X11 = load([directory,'x_p1_c1.txt']);
X12 = load([directory,'x_p1_c2.txt']);
X14 = load([directory,'x_p1_c4.txt']);
X18 = load([directory,'x_p1_c8.txt']);

%Get data for p = 2
U21 = load([directory,'u_p2_c1.txt']);
U22 = load([directory,'u_p2_c2.txt']);
U24 = load([directory,'u_p2_c4.txt']);
U28 = load([directory,'u_p2_c8.txt']);
Q21 = load([directory,'q_p2_c1.txt']);
Q22 = load([directory,'q_p2_c2.txt']);
Q24 = load([directory,'q_p2_c4.txt']);
Q28 = load([directory,'q_p2_c8.txt']);
X21 = load([directory,'x_p2_c1.txt']);
X22 = load([directory,'x_p2_c2.txt']);
X24 = load([directory,'x_p2_c4.txt']);
X28 = load([directory,'x_p2_c8.txt']);

%Get data for p = 3
U31 = load([directory,'u_p3_c1.txt']);
U32 = load([directory,'u_p3_c2.txt']);
U34 = load([directory,'u_p3_c4.txt']);
U38 = load([directory,'u_p3_c8.txt']);
Q31 = load([directory,'q_p3_c1.txt']);
Q32 = load([directory,'q_p3_c2.txt']);
Q34 = load([directory,'q_p3_c4.txt']);
Q38 = load([directory,'q_p3_c8.txt']);
X31 = load([directory,'x_p3_c1.txt']);
X32 = load([directory,'x_p3_c2.txt']);
X34 = load([directory,'x_p3_c4.txt']);
X38 = load([directory,'x_p3_c8.txt']);


%% Get the figures

    %For p = 0
    figure
    plot(X01,U01(2:end),X02,U02(2:end),'rx-',X04,U04(2:end),'g.-',X08,U08(2:end))
    xlabel('Position (m)')
    ylabel('Concentration (kg/m^3)')
    title('Concentration vs Position at t = 2 seconds for p = 0')
    legend('Nx = 10','Nx = 20','Nx = 40','Nx = 80',0)
    axis([0,2*pi,1.1*min(U08(2:end)),1.1*max(U08(2:end))])
    
    figure
    plot(X01,Q01(2:end),X02,Q02(2:end),'rx-',X04,Q04(2:end),'g.-',X08,Q08(2:end))
    xlabel('Position (m)')
    ylabel('Diffusive Flux (kg/(m^2s))')
    title('Diffusive Flux vs Position at t = 2 seconds for p = 0')
    legend('Nx = 10','Nx = 20','Nx = 40','Nx = 80',0)
    axis([0,2*pi,1.1*min(Q08(2:end)),1.1*max(Q08(2:end))])
    
    %For p = 1
    figure
    plot(X11,U11(2:end),X12,U12(2:end),'rx-',X14,U14(2:end),'g.-',X18,U18(2:end))
    xlabel('Position (m)')
    ylabel('Concentration (kg/m^3)')
    title('Concentration vs Position at t = 2 seconds for p = 1')
    legend('Nx = 10','Nx = 20','Nx = 40','Nx = 80',0)
    axis([0,2*pi,1.1*min(U18(2:end)),1.1*max(U18(2:end))])
    
    figure
    plot(X11,Q11(2:end),X12,Q12(2:end),'rx-',X14,Q14(2:end),'g.-',X18,Q18(2:end))
    xlabel('Position (m)')
    ylabel('Diffusive Flux (kg/(m^2s))')
    title('Diffusive Flux vs Position at t = 2 seconds for p = 1')
    legend('Nx = 10','Nx = 20','Nx = 40','Nx = 80',0)
    axis([0,2*pi,1.1*min(Q18(2:end)),1.1*max(Q18(2:end))])
    
        %For p = 2
    figure
    plot(X21,U21(2:end),X22,U22(2:end),'rx-',X24,U24(2:end),'g.-',X28,U28(2:end))
    xlabel('Position (m)')
    ylabel('Concentration (kg/m^3)')
    title('Concentration vs Position at t = 2 seconds for p = 2')
    legend('Nx = 10','Nx = 20','Nx = 40','Nx = 80',0)
    axis([0,2*pi,1.1*min(U28(2:end)),1.1*max(U28(2:end))])
    
    figure
    plot(X21,Q21(2:end),X22,Q22(2:end),'rx-',X24,Q24(2:end),'g.-',X28,Q28(2:end))
    xlabel('Position (m)')
    ylabel('Diffusive Flux (kg/(m^2s))')
    title('Diffusive Flux vs Position at t = 2 seconds for p = 2')
    legend('Nx = 10','Nx = 20','Nx = 40','Nx = 80',0)
    axis([0,2*pi,min(Q28(2:end))*1.1,1.1*max(Q28(2:end))])
    
    %For p = 3
    figure
    plot(X31,U31(2:end),X32,U32(2:end),'rx-',X34,U34(2:end),'g.-',X38,U38(2:end))
    xlabel('Position (m)')
    ylabel('Concentration (kg/m^3)')
    title('Concentration vs Position at t = 2 seconds for p = 3')
    legend('Nx = 10','Nx = 20','Nx = 40','Nx = 80',0)
    axis([0,2*pi,1.1*min(U38(2:end)),1.1*max(U38(2:end))])
    
    figure
    plot(X31,Q31(2:end),X32,Q32(2:end),'rx-',X34,Q34(2:end),'g.-',X38,Q38(2:end))
    xlabel('Position (m)')
    ylabel('Diffusive Flux (kg/(m^2s))')
    title('Diffusive Flux vs Position at t = 2 seconds for p = 3')
    legend('Nx = 10','Nx = 20','Nx = 40','Nx = 80',0)
    axis([0,2*pi,1.1*min(Q38(2:end)),1.1*max(Q38(2:end))])
