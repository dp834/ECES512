%% Question 2
% System and pole definition

A = [ 1 1 -2 ;
      0 1  1 ;
      0 0  1];
B = [ 1 ;
      0 ;
      1];
poles = [-2 -1+i -1-i];

% Using place
gains_place = place(A,B,poles);
gains_place
% gains_place =
% 15.0000   47.0000   -8.0000



% Using Acker
gains_acker = acker(A,B,poles);
gains_acker

% gains_acker =
% 15    47    -8


% lsim the systems
sys_orig = ss(A,B,[1 0 0],0);
sys_modi = ss(A-B*gains_place,B,[1 0 0],0);
t = 0:.01:10;
u = t >= 0;

[~,~,xs] = lsim(sys_orig, u,t);
plot(t,xs);
legend("x1","x2","x3");
saveas(gcf, "images/p2_original_system.png");

[~,~,xs] = lsim(sys_modi, u,t);
plot(t,xs);
legend("x1","x2","x3");
saveas(gcf, "images/p2_modified_system.png");
