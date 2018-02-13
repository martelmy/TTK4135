load('pitch.mat');
load('travel.mat');
load('travel_ref.mat');
load('pitch_ref.mat');

figure(1); clf
hold on;
grid on;
xlabel('time [s]');
ylabel('pitch [rad]');
y = measured_pitch(2,:);
t = measured_pitch(1,:);
p1 = plot(t,y, 'b');
y2 = pitch_ref(2,:);
t2 = pitch_ref(1,:);
p2 = plot(t2,y2,'r');
title('Pitch trajectory with Q = diag([100 0.1 10 0.1])');
legend([p1,p2], 'pitch', 'pitch reference');
hold off;
%saveas(gcf, '10_3_2pitchQ3','epsc');

figure(2); clf
hold on;
grid on;
xlabel('time [s]');
ylabel('travel [rad]');
y1 = measured_travel(2,:);
t1 = measured_travel(1,:);
p3 = plot(t1,y1,'b');
y3 = travel_ref(2,:);
t3 = travel_ref(1,:);
p4 = plot(t3,y3,'r');
title('Travel trajectory with Q = diag([100 0.1 10 0.1])');
legend([p3,p4], 'travel', 'travel reference');
hold off;

%saveas(gcf, '10_3_2travelQ3','epsc');