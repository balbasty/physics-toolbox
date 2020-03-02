%% EXP MODEL

obs = 0.5;
xhat = log(obs);
f = @(x) (exp(x)-obs).^2;  % Objective function
g = @(x) 2*(exp(x)-obs).*exp(x); % Gradient
H = @(x) 2*exp(2*x);  % Approximate Hessian
b = @(x0,x) f(x0) + g(x0).*(x-x0) + 0.5*H(x0).*(x-x0).^2;

x = linspace(-10,5,1024);
figure
plot(x,f(x));
lim = minmax(f(x));
hold on
plot(x,b(xhat,x));
plot(x,b(-5,x));
plot(x,b(0,x));
hold off
ylim([0 0.5]);

%% EXPoEXP MODEL

obs = 0.5;
xhat = log(-log(obs));
f = @(x) (exp(-exp(x))-obs).^2;  % Objective function
g = @(x) 2*(exp(-exp(x))-obs).*(-exp(x-exp(x))); % Gradient
H = @(x) 2*(exp(2*x-2*exp(x)));  % Approximate Hessian
b = @(x0,x) f(x0) + g(x0).*(x-x0) + 0.5*H(x0).*(x-x0).^2;
H2 = @(x) 2*(exp(2*x));  % Approximate Hessian
b2 = @(x0,x) f(x0) + g(x0).*(x-x0) + 0.5*H2(x0).*(x-x0).^2;

x = linspace(-10,5,1024);
figure
plot(x,f(x));
lim = minmax(f(x));
hold on
plot(x,b(xhat,x));
plot(x,b(-5,x));
plot(x,b(5,x));
plot(x,b2(xhat,x));
plot(x,b2(-5,x));
plot(x,b2(2,x));
hold off
ylim([0 0.5]);