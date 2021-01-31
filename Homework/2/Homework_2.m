%% Constants
global R L m K g
R = 1;
L = 0.01;
m = 0.05;
K = 0.0001;
g = 9.81;

%% Script
% Finding equilibrium points
for i = [7]
    fprintf("Equilibrium for u(t) = %f\n", i)
    tmp = equilibrium_constant(i);
    fprintf("x_1: %g + %gi\nx_2: %g + %gi\nx_3: %g + %gi\n", real(tmp(1)), imag(tmp(1)),real(tmp(2)), imag(tmp(2)),real(tmp(3)), imag(tmp(3)));
end


% For accessing the linear system cells
linear_systems = {};
for i = [7]
    fprintf("Linearization for u(t) = %f\n", i)
    [A, B, C, D] = linearize(i);
    linear_systems{end+1} = {ss(A, B, C, D), equilibrium_constant(i), i};
end

t = 0:0.0005:.2;
u = zeros(size(t));
system_conditions = {{7, [7 -.010 0]}};
for cond = system_conditions
    for sys = linear_systems
        % convert coordinates to the local linear approx coordintaes u(:,:) = cond{1}{1} - sys{1}{3};
		u(:,:) = cond{1}{1}-sys{1}{3};
        [~, ~, y] = lsim(sys{1}{1}, u, t, cond{1}{2} - sys{1}{2}.');
        [~, yode] = ode45(@(t,x) dx(x,@(x) cond{1}{1}, t), t, cond{1}{2});
        y = y + sys{1}{2}.';
        subplot(3,1,1);
        plot(t,y(:,1));
		hold on;
        plot(t,yode(:,1),"k");
		fplot(7, [t(1) t(end)], "b*");
		hold off;
        xlabel("Time");
        ylabel("Current in Amps");
        subplot(3,1,2);
        plot(t,y(:,2));
		hold on;
        plot(t,yode(:,2),"k");
		fplot(@(t) (0.199807E-1+(-0.499033E-2)*exp((-0.313369E2)*t)+(-0.499033E-2)*exp(0.313369E2*t)), [t(1) t(end)], "b*");
		hold off;
        xlabel("Time");
        ylabel("Distance in meters");
        subplot(3,1,3);
        plot(t,y(:,3));
		hold on;
        plot(t,yode(:,3),"k");
		fplot(@(t) (0.156381E0*exp((-0.313369E2)*t)+(-0.156381E0)*exp(0.313369E2*t)), [t(1) t(end)], "b*");
		hold off;
        xlabel("Time");
        ylabel("Velocity in m/s");
        plot_title = sprintf("lsim of system linearized around input = %f \n u(t) = %f  x_0 = %f", sys{1}{3}, cond{1}{1}, cond{1}{2}(2));
        sgtitle(plot_title);
        filename = sprintf("images/linearized_%f_input_%f_initial_pos_%f.png", sys{1}{3}, cond{1}{1}, cond{1}{2}(2));
        saveas(gcf, filename);
    end
end
clf();


%% Functions
function x_dot = dx(x,func_impulse, t)
    global R L m K g
    % dx_1 = (u(t) - R*x_1(t))/(L);
    % dx_2 = x_3(t);
    % dx_3 = m*g - (K*x_1^2(t))/(x_2(t));
    x_dot = [
        (func_impulse(t) - R*x(1))/(L);
        x(3);
        m*g - (K*x(1)^2)/(x(2));
    ];
end

% returns the equilibrium given a constant input function
function eq = equilibrium_constant(u_t)
    global R L m K g
    eq = [
    u_t/R;
    K*u_t^2/(m*g*R^2);
    0;
    ];
end

% Linearize given an constant input function
function [A, B, C, D] = linearize(u_t)
    global R L m K g
    eq = equilibrium_constant(u_t);
    A = [
    -R/L 0 0;
    0 0 1;
    -2*K*eq(1)/(m*eq(2)) (K/m)*(eq(1)/(eq(2)))^2 0;
    ];
    B = [
    1/L;
    0;
    0;
    ];
    C = [
    0 0 0;
    0 1 0;
    0 0 0;
    ];
    D = [0];
end
