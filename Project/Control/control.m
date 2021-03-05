%% Constants
global R L m K g
R = 1;
L = 0.01;
m = 0.05;
K = 0.0001;
g = 9.81;


%% Full State Feedback
A = [ -100  0  0;
        0   0  1;
    -2.803 982 0];
B = [1/.01 ; 0 ; 0];
C = [0 1 0];
D = 0;
x0 = [0 0 0];
linearization_point = [7 .00998 0];

t = 0:.001:.5;

sys_orig = ss(A, B, C, D);
u = 0* (t>0);
[~,~,xs] = lsim(sys_orig, u,t, x0 - linearization_point);
xs = xs + linearization_point;
subplot(3,1,1);
sgtitle(" FullState Feedback Reference = 0");
plot(t,xs(:,1));
legend("Voltage");
subplot(3,1,2);
plot(t,xs(:,2));
legend("Distance")
subplot(3,1,3);
plot(t,xs(:,3));
legend("Velocity")
saveas(gcf, "images/OpenLoop_No input.png");

kPoles = [ -20+20i -20-20i -100];

ctrb(A,B)
k = acker(A,B,kPoles)

sys_fsf = ss(A-B*k, B, C, D);

u = 0* (t>0);
[~,~,xs] = lsim(sys_fsf, u,t, x0 - linearization_point);
xs = xs + linearization_point;
subplot(3,1,1);
sgtitle(" FullState Feedback Reference = 0");
plot(t,xs(:,1));
legend("Voltage");
subplot(3,1,2);
plot(t,xs(:,2));
legend("Distance")
subplot(3,1,3);
plot(t,xs(:,3));
legend("Velocity")
saveas(gcf, "images/FullState_Feedback_ref_0.png");

u = 0.005* (t>0);
[~,~,xs] = lsim(sys_fsf, u,t, x0 - linearization_point);
xs = xs + linearization_point;
subplot(3,1,1);
sgtitle(" FullState Feedback with no NBar correction, Reference = .015");
plot(t,xs(:,1));
legend("Voltage");
subplot(3,1,2);
plot(t,xs(:,2));
legend("Distance")
subplot(3,1,3);
plot(t,xs(:,3));
legend("Velocity")
saveas(gcf, "images/FullState_Feedback_No_nBar_ref_p015.png");

nbar = -285.72;

u = 0.005*(t>0);
[~,~,xs] = lsim(sys_fsf, nbar*u,t, x0 - linearization_point);
xs = xs + linearization_point;
subplot(3,1,1);
sgtitle(" FullState Feedback with NBar correction, Reference = .015");
plot(t,xs(:,1));
legend("Voltage");
subplot(3,1,2);
plot(t,xs(:,2));
legend("Distance")
subplot(3,1,3);
plot(t,xs(:,3));
legend("Velocity")
saveas(gcf, "images/FullState_Feedback_with_nBar_ref_p015.png");

lPoles = [-100+100i -100-100i -500];
obsv(A,C)
l = place(A',C', lPoles)
l = l';
A_of = [A -B*k;
        l*C A-B*k-l*C];
B_of = [B;B];
C_of = [C 0 0 0; 0 0 0 C];
D_of = 0;

sys_of = ss(A_of, B_of, C_of, D_of);

u = 0* (t>0);

[~,~,xs] = lsim(sys_of, nbar*u,t, [(x0-linearization_point) 0 0 0]);
xs(:,1:3) = xs(:,1:3) + linearization_point;
xs(:,4:6) = xs(:,4:6) + linearization_point;

subplot(3,1,1);
sgtitle("Output Feedback Reference = 0");
plot(t,xs(:,1));
hold on;
plot(t,xs(:,4),"--");
plot(t,xs(:,4)-xs(:,1));
hold off
title("Voltage")
legend("Actual", "Estimated", "Error");
subplot(3,1,2);
plot(t,xs(:,2));
hold on;
plot(t,xs(:,5),"--");
plot(t,xs(:,5)-xs(:,2));
hold off
title("Distance")
legend("Actual", "Estimated", "Error");
subplot(3,1,3);
plot(t,xs(:,3));
hold on;
plot(t,xs(:,6),"--");
plot(t,xs(:,6)-xs(:,3));
hold off
title("Velocity")
legend("Actual", "Estimated", "Error");
saveas(gcf, "images/Output_Feedback_ref_0.png");

u = 0.005* (t>0);

[~,~,xs] = lsim(sys_of, nbar*u,t, [(x0-linearization_point) 0 0 0]);
xs(:,1:3) = xs(:,1:3) + linearization_point;
xs(:,4:6) = xs(:,4:6) + linearization_point;

subplot(3,1,1);
sgtitle("Output Feedback Reference = 0.015");
plot(t,xs(:,1));
hold on;
plot(t,xs(:,4),"--");
plot(t,xs(:,4)-xs(:,1));
hold off
title("Voltage")
legend("Actual", "Estimated", "Error");
subplot(3,1,2);
plot(t,xs(:,2));
hold on;
plot(t,xs(:,5),"--");
plot(t,xs(:,5)-xs(:,2));
hold off
title("Distance")
legend("Actual", "Estimated", "Error");
subplot(3,1,3);
plot(t,xs(:,3));
hold on;
plot(t,xs(:,6),"--");
plot(t,xs(:,6)-xs(:,3));
hold off
title("Velocity")
legend("Actual", "Estimated", "Error");
saveas(gcf, "images/Output_Feedback_ref_p015.png");



%% Messing around with ode45 to see how the linearized system works
% [ts, y] = ode45(@(t,x) dx(x, @(x) cond{1}{1}, t), [t(1) t(end)], cond{1}{2});

%% Functions
function x_dot = step(x, t, linear_point, A, B, C, K, L)
    force = -B(3:end,1:3)*K*x(3:6);
end

function x_dot = dx_nonLinear(x,input)
    global R L m K g
    % dx_1 = (u(t) - R*x_1(t))/(L);
    % dx_2 = x_3(t);
    % dx_3 = m*g - (K*x_1^2(t))/(x_2(t));

    x_dot = [
        (input - R*x(1))/(L);
        x(3);
        m*g - (K*x(1)^2)/(x(2));
    ];
end
