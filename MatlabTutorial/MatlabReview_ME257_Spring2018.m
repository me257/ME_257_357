clear; close all; clc
%% Basics
% how are numbers represented in Matlab: vectors, matrices and their
% operations

% useful functions that help you to creat a vector: zeros, ones, linspace
a = zeros(2, 3);
b = ones(2, 3);
c = linspace(0, 1, 11);

% creating a struct
% e.g. a struct named pizza that has the following components:
% number : the index of the pizza
% deliver : time each pizza is delivered
% receive : time each pizza is received
% method 1:
pizza = struct('number', linspace(0, 5, 6), ...
    'deliver', zeros(1, 6), ...
    'receive', ones(1, 6));

% method 2:
pizza.number = linspace(0, 5, 6);
pizza.deliver = zeros(1, 6);
pizza.receive = zeros(1, 6);

%% Control statements
% loops: for, while
for i = 1 : 6
    pizza.deliver(i) = rand(1) + pizza.number(i);
    fprintf('Peter delievers the %i th pizza at %5.3f\n', ...
        pizza.number(i), pizza.deliver(i));
end

i = 1;
while i <= 6
    pizza.receive(i) = rand(1) + pizza.deliver(i);
    fprintf('Hao receives the %i th pizza at %5.3f\n', i, pizza.receive(i));
    i = i + 1;
end

% note:
% operations above can be performed using the pointwise operations
% (equivalent but faster)
% pizza.deliver = rand(1, 6) + pizza.number;
% pizza.receive = rand(1, 6) + pizza.deliver;

% conditional statements: if, switch
i = '10';
if (strcmp(i, '8'))
    fprintf('hello!\n');
elseif (strcmp(i, '9'))
    fprintf('hi\n');
else
    fprintf('bye\n');
end
    
switch i
    case '8'
        fprintf('hello!\n');
    case '9'
        fprintf('hi\n');
    otherwise
        fprintf('bye\n');
end
%% Functions
% write function in a seperate file

% write a lambda function
dydt = @(t, x) [x(2); x(1) * exp(-0.5 * t)];

%% Solve ODE (systems)
% Problem:
% dx / dt = u
% du / dt = a
% where a = x * exp(-0.5 * t)

res = ode45(dydt, [0 20], [0; 1]);
% res = ode23s(dydt, [0 20], [0; 1]);

%% Plots
figure
l = plot(res.x, res.y(2, :));
set(gca, 'FontSize', 16)
set(l, 'LineWidth', 4, 'color', 'k')
xlabel('$t$ [s]', 'Interpreter','latex','FontSize',16)
ylabel('$u$ [m/s]', 'Interpreter','latex','FontSize',16)
leg = legend('velocity');
set(leg, 'Interpreter','latex', 'FontSize',16, 'location','NorthWest')