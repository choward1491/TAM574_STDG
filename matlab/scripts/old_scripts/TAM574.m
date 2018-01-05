function pdot = TAM574(t,x,Dp,Dm,kappa,tau,a)

c = sqrt(kappa/tau);
l1 = a+c;
l2 = a-c;
len = length(x);
u = x(1:len/2);
q = x(len/2+1:end);



Ap = .5.*[l1,l1/c;c*l1,l1];
Am = .5.*[l2,-l2/c;-c*l2,l2];
D = [0,0;0,1/tau];

dudxp = Dp*u; 
dudxm = Dm*u;
dqdxp = Dp*q;
dqdxm = Dm*q;

qdot = zeros(len/2,1);
udot = qdot;

for i = 1:len/2
    
    tdot = -(Ap*[dudxp(i);dqdxp(i)] + Am*[dudxm(i);dqdxm(i)] + D*[u(i);q(i)]);
    
    udot(i) = tdot(1);
    qdot(i) = tdot(2);
    
end

pdot = [udot;qdot];
end