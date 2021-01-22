%% Problem 1
A = [-2 0; 1 -1];
b = [1;0];
ts = (0:.01:10)';
u = ts>=0;
x0 = [2;3];
sys = ss(A,b,[1 0],0);

[~,ts,x] = lsim(sys, u, ts, x0);

plot(ts, x);
hold on;
fplot(@(x) ( (3/2)*exp(-2*x)+ (1/2)), [ts(1) ts(end)], "k*");
fplot(@(x) ( (-3/2)*exp(-2*x) + 4*exp(-x) + (1/2)), [ts(1) ts(end)], "g*");
legend("x1 Sim", "x2 Sim", "x1 Prediction", "x2 Prediction");
saveas(gcf, "images/p1_predicted_vs_simulation.png");

hold off;
cla

%% Problem 2
sys = tf([3 2],[1 3 2 0]);
ts = (0:.01:10)';
u = ts>=0;
sys
[y, ts] = lsim(sys, u, ts);
plot(ts, y);
hold on;
fplot(@(x)  (exp(-2*x)-exp(-x)+x), [ts(1) ts(end)], "b*");
legend("Sim", "Prediction");
saveas(gcf, "images/p2_2_predicted_vs_simulation.png");

hold off;
cla
