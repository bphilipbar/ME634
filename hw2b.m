
stopEta = 5;
pointsEta = 51;
iterations = 20;



h = stopEta / (pointsEta - 1);  % step size

eta = linspace(0,stopEta,pointsEta)';

% define our starting guesses (see Mathematica file for plots)
% you can verify the boundary conditions are true here
f1 = eta;
f2 = 1 - exp(-eta);
f3 = (1-tanh(eta-3))/6;



for j = 1:iterations
    
    % evaluate our guesses/solutions at their centered values
    halfF1 = (f1(2:end) + f1(1:end-1))/2;
    halfF2 = (f2(2:end) + f2(1:end-1))/2;
    halfF3 = (f3(2:end) + f3(1:end-1))/2;
    
    % loop over all points with equations to build up the matrix and b vector
    % (see Mathematica file for the equations that gave us these)
    matrix = [];
    b = [];
    for i = 1:(pointsEta - 1)
        A = [-1/h, -1/2, 0; 0, -1/h, -1/2; halfF3(i)/4, 0 , halfF1(i)/4-1/h];
        B = [1/h, -1/2, 0; 0, 1/h, -1/2; halfF3(i)/4, 0 , halfF1(i)/4+1/h];
        row = [repmat(zeros(3),1,i-1),A,B,repmat(zeros(3),1,pointsEta-1-i)];
        matrix = [matrix;row];
        b=[b;-(f1(i+1)-f1(i))/h+halfF2(i);-(f2(i+1)-f2(i))/h+halfF3(i); ...
            -(f3(i+1)-f3(i))/h-halfF1(i)*halfF3(i)/2];
    end
    
    % add the final three rows to give boundary conditions
    row1 = zeros(1,pointsEta*3);
    row2 = zeros(1,pointsEta*3);
    row3 = zeros(1,pointsEta*3);
    row1(1) = 1;  % df1(1) = 0
    row2(2) = 1;  % df2(1) = 0
    row3(end-1) = 1;  % df2(end) = 0
    b = [b;0;0;0];
    matrix = [matrix;row1;row2;row3];
    
    % solve the system
    sol = matrix\b;
    df1 = sol(1:3:pointsEta*3);
    df2 = sol(2:3:pointsEta*3);
    df3 = sol(3:3:pointsEta*3);
    
    % get new iteration
    f1 = f1 + df1;
    f2 = f2 + df2;
    f3 = f3 + df3;
    
end

% sometimes round-off creeps in to give a negative value screwing up the
% plots...
f1(1) = 0;
f2(1) = 0;

% make plots
close all
figure('position', [0, 0, 1800, 700])
subplot(1,3,1)
plot(eta,f1,'r')
xlabel('eta'), ylabel('f')
subplot(1,3,2)
plot(eta,f2,'r')
xlabel('eta'), ylabel('f prime')
subplot(1,3,3)
plot(eta,f3,'r')
xlabel('eta'), ylabel('f double prime')

