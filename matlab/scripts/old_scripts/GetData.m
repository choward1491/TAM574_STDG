%Read in data

directory = '../data/';
U0 = load([directory,'u_p1_nx100_c256.txt']);
Q0 = load([directory,'q_p1_nx100_c256.txt']);

x = load([directory,'xvector.txt']);
x_nx100 = load([directory,'x_p1_nx100_c256.txt']);

t = U0(:,1);
u0 = U0(:,2:end); clear U0;
q0 = Q0(:,2:end); clear Q0;

sprintf('I am done yo')

%% Plot some things
figure
for i = 1:length(t)
    subplot(2,1,1)
    plot(x,u0(i,:),x,u3(i,:))
    
    subplot(2,1,2)
    plot(x,q0(i,:),x,q3(i,:))
    
    pause(0.01)
end

%% Plot some independent things
uexact = @(x,t) exp(-t).*sin(x);
index = 1:length(t);
tw = 1:5;

for i = index(t == 2)
    
    figure
    subplot(2,1,1)
    plot(x_nx100,u0(i,:),x_nx100,uexact(x_nx100,t(i)))
    
    subplot(2,1,2)
    plot(x_nx100,u0(i,:))
    pause(0.01)
end

%% Get just the values desired for plotting
U0 = load([directory,'u_p1_nx100_c256.txt']);
Q0 = load([directory,'q_p1_nx100_c256.txt']);
t = U0(:,1);
u0 = U0(:,2:end); clear U0;
q0 = Q0(:,2:end); clear Q0;
x = load([directory,'x_p1_nx100_c256.txt']);
index = 1:length(t);

Uw = zeros(5,length(x));
Qw = Uw;

for i = 1:5
Uw(i,:) = u0(index( t == i ), :);
Qw(i,:) = q0( index( t == i ),:);
end

save('u_nx100.txt','Uw','-ASCII');
save('q_nx100.txt','Qw','-ASCII');
save('x_nx100.txt','x','-ASCII')

clear
end_msg = sprintf('The work is done my Commander')