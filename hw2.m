%% have MATLAB solve the equation

guess = 0.332;
finalEta = 5; %should be large when modifying the above value

f = @(t,x) [x(2);x(3);-x(1)*x(3)/2];
[etaa,soll] = ode45(f,[0,finalEta],[0,0,guess]);

soll(end,2) %check final f' value


%% plots

close all
plot(etaa,soll(:,1))
xlabel('eta'), ylabel('f')
figure
plot(etaa,soll(:,2))
xlabel('eta'), ylabel('f prime')
figure
plot(etaa,soll(:,3))
xlabel('eta'), ylabel('f double prime')


%% plots for comparison if hw2b.m has already been run

% red is the Keller Box Method (black is MATLAB's solution)
close all
figure('position', [0, 0, 1800, 700])
subplot(1,3,1)
plot(etaa,soll(:,1),'k',eta,f1,'r')
xlabel('eta'), ylabel('f')
subplot(1,3,2)
plot(etaa,soll(:,2),'k',eta,f2,'r')
xlabel('eta'), ylabel('f prime')
subplot(1,3,3)
plot(etaa,soll(:,3),'k',eta,f3,'r')
xlabel('eta'), ylabel('f double prime')
