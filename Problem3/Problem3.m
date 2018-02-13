init;
problem2;
%Q_lq = diag([100 10 10 0.1]);
Q_lq = diag([4 4 1 1]); %optimal
R_lq = diag(1);

[K, S, e] = dlqr(A1,B1,Q_lq,R_lq);

t1 = [1:delta_t:36]';
x_t = [t1 x1 x2 x3 x4];
u_t = [t1 u];
x1_t = [t1 x1];