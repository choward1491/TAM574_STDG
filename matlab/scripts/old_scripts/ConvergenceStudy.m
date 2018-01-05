%Convergence study stuff
close all
Nx = [10,20,40,80];
h = 2*pi./Nx;

Zeta_U = [2.19641,3.19704,4.21895,5.0508;
               1.23482,0.721079,0.403697,0.215252;
               0.1427,0.0374825,0.00977865,0.00251122;
                0.00819868,0.00108047,0.000139279,1.77078E-05];
            
Zeta_Q = [55.6964,80.7534,108.344,130.583;
               10.7452,8.40573,5.55747,3.21516;
               1.74891,0.498069,0.140258,0.037682;   
                0.0987838,0.0143293,0.00196368,0.000258194];

%% Find the log of the zetas and h

LZU = log(Zeta_U)/log(10);
LZQ = log(Zeta_Q)/log(10);

Lh = log(h)/log(10);            
            
%% Plot the figures log log
figure 
loglog(h,Zeta_U(1,:),h,Zeta_U(2,:), h,Zeta_U(3,:), h, Zeta_U(4,:) )
title('\zeta^u vs Element Spacial Length')
xlabel('Element Spacial Length')
ylabel('\zeta^u')
legend('P = 0','P = 1','P = 2','P = 3',0)

figure 
loglog(h,Zeta_Q(1,:),h,Zeta_Q(2,:), h,Zeta_Q(3,:), h, Zeta_Q(4,:) )
title('\zeta^q vs Element Spacial Length')
xlabel('Element Spacial Length')
ylabel('\zeta^q')
legend('P = 0','P = 1','P = 2','P = 3',0)

%% Plot the figures normally
figure 
plot(Lh,LZU(1,:),Lh,LZU(2,:), Lh,LZU(3,:), Lh, LZU(4,:) )
title('log(\zeta^u) vs log(h)')
xlabel('log(h)')
ylabel('log(\zeta^u)')
legend('P = 0','P = 1','P = 2','P = 3',0)

figure 
plot(Lh,LZQ(1,:),Lh,LZQ(2,:), Lh,LZQ(3,:), Lh, LZQ(4,:) )
title('log(\zeta^q) vs log(h)')
xlabel('log(h)')
ylabel('log(\zeta^q)')
legend('P = 0','P = 1','P = 2','P = 3',0)
%% Find the convergence rates
p0u = ( LZU(1,end)-LZU(1,1) )/ ( Lh(1,end) - Lh(1,1) );
p0q = ( LZQ(1,end)-LZQ(1,1) )/ ( Lh(1,end) - Lh(1,1) );

p1u = ( LZU(2,end)-LZU(2,1) )/ ( Lh(1,end) - Lh(1,1) );
p1q = ( LZQ(2,end)-LZQ(2,1) )/ ( Lh(1,end) - Lh(1,1) );

p2u = ( LZU(3,end)-LZU(3,1) )/ ( Lh(1,end) - Lh(1,1) );
p2q = ( LZQ(3,end)-LZQ(3,1) )/ ( Lh(1,end) - Lh(1,1) );

p3u = ( LZU(4,end)-LZU(4,1) )/ ( Lh(1,end) - Lh(1,1) );
p3q = ( LZQ(4,end)-LZQ(4,1) )/ ( Lh(1,end) - Lh(1,1) );

Convergence = [p0u,p0q;p1u,p1q;p2u,p2q;p3u,p3q]