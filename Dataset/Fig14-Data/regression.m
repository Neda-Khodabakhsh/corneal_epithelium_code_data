% داده‌های فرضی
t = 0:5:1500; % زمان
d1 = 500 - matrix_tip001'; % مسافت برای r_shed = 0.001
d2 = 500 - matrix_tip005'; % مسافت برای r_shed = 0.005
d3 = 500 - matrix_tip009'; % مسافت برای r_shed = 0.009
d4 = 500 - matrix_tip01'; % مسافت برای r_shed = 0.01
d5 = 500 - matrix_tip015'; % مسافت برای r_shed = 0.015
d6 = 500 - matrix_tip02'; % مسافت برای r_shed = 0.02

% انجام رگرسیون خطی
p1 = polyfit(t, d1, 1); % برازش خطی برای r_shed = 0.001
p2 = polyfit(t, d2, 1); % برازش خطی برای r_shed = 0.005
p3 = polyfit(t, d3, 1); % برازش خطی برای r_shed = 0.009
p4 = polyfit(t, d4, 1); % برازش خطی برای r_shed = 0.01
p5 = polyfit(t, d5, 1); % برازش خطی برای r_shed = 0.015
p6 = polyfit(t, d6, 1); % برازش خطی برای r_shed = 0.02

% استخراج سرعت میانگین (شیب خط)
v1 = p1(1); 
v2 = p2(1);
v3 = p3(1);
v4 = p4(1);
v5 = p5(1);
v6 = p6(1);


% نمایش نتایج
fprintf('Average velocity for r_shed = 0.001: %.2f\n', v1);
fprintf('Average velocity for r_shed = 0.005: %.2f\n', v2);
fprintf('Average velocity for r_shed = 0.009: %.2f\n', v3);
fprintf('Average velocity for r_shed = 0.01: %.2f\n',  v4);
fprintf('Average velocity for r_shed = 0.015: %.2f\n', v5);
fprintf('Average velocity for r_shed = 0.02: %.2f\n',  v6);

% رسم نمودار
figure;
hold on;
scatter(t, d1, 'bo', 'DisplayName', 'Data (r\_shed = 0.001)');
scatter(t, d2, 'ro', 'DisplayName', 'Data (r\_shed = 0.005)');
scatter(t, d3, 'ro', 'DisplayName', 'Data (r\_shed = 0.009)');
scatter(t, d4, 'ro', 'DisplayName', 'Data (r\_shed = 0.01)');
scatter(t, d5, 'ro', 'DisplayName', 'Data (r\_shed = 0.015)');
scatter(t, d6, 'ro', 'DisplayName', 'Data (r\_shed = 0.02)');



plot(t, polyval(p1, t), 'b--', 'DisplayName', sprintf('Fit (v1 = %.2f)', v1));
plot(t, polyval(p2, t), 'r--', 'DisplayName', sprintf('Fit (v2 = %.2f)', v2));
plot(t, polyval(p3, t), 'r--', 'DisplayName', sprintf('Fit (v3 = %.2f)', v3));
plot(t, polyval(p4, t), 'r--', 'DisplayName', sprintf('Fit (v4 = %.2f)', v4));
plot(t, polyval(p5, t), 'r--', 'DisplayName', sprintf('Fit (v5 = %.2f)', v5));
plot(t, polyval(p6, t), 'r--', 'DisplayName', sprintf('Fit (v6 = %.2f)', v6));

xlabel('Time (t)');
ylabel('Distance (d)');
title('Linear Regression Fit for Distance vs. Time');
legend;
grid on;
hold off;
