%Get all the data for plotting and make the plots at times
% t = 1,2,3,4,5

close all
clear all

directory = '../data/';
%Get data for C = 1/2
Uc2 = load([directory,'u_c2.txt']);
%Qc2 = load([directory,'q_c2.txt']);
Xc2 = load([directory,'x_c2.txt']);

%Get data for C = 1/4
Uc4 = load([directory,'u_c4.txt']);
%Qc4 = load([directory,'q_c4.txt']);
Xc4 = load([directory,'x_c4.txt']);

%Get data for C = 1/16
Uc16 = load([directory,'u_c16.txt']);
%Qc16 = load([directory,'q_c16.txt']);
Xc16 = load([directory,'x_c16.txt']);

%Get data for C = 1/256
Uc256 = load([directory,'u_c256.txt']);
%Qc256 = load([directory,'q_c256.txt']);
Xc256 = load([directory,'x_c256.txt']);

%Get data for C = 1/256 and Nx = 100
Unx = load([directory,'u_nx100.txt']);
%Qnx = load([directory,'q_nx100.txt']);
Xnx = load([directory,'x_nx100.txt']);

uexact = @(x,t) exp(-t).*sin(x);
x = linspace(0,2*pi,200);

t = 1:5;

for i = 1:length(t)
    time_str = num2str(t(i));
    U = uexact(x,t(i));
   minu = min(U);
   maxu = max(U);
    
    figure
    plot(Xc2,Uc2(i,:),Xc4,Uc4(i,:),Xc16,Uc16(i,:),Xc256,Uc256(i,:),Xnx,Unx(i,:),x,U)
    title(['Concentration vs Position at t = ',time_str,' seconds for p = 1'])
    axis([0,2*pi,minu*1.1,maxu*1.1])
    xlabel('Position (m)')
    ylabel('Concentration (kg/m^3)')
    legend('C = 2^{-1}','C = 2^{-2}','C = 4^{-2}','C = 16^{-2}','C = 16^{-2} & Nx = 100','Exact Parabolic Soln',0)
    
end