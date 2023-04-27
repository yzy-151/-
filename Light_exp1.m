
%%
%函数的基本公式与参数
R = 1;
imax = 0.5;
k1 = 2;k2 = 3;k3 = 50;
vled = -0.5:0.1:1;

%矩阵要用.^,除法的时候也是，当下面有除数时要用点
%求公式中的h
f_vled = vled/R;
h_vled1 = f_vled./((1+(f_vled/imax).^(2*k1)).^(0.5/k1));
h_vled2 = f_vled./((1+(f_vled/imax).^(2*k2)).^(0.5/k2));
h_vled3 = f_vled./((1+(f_vled/imax).^(2*k3)).^(0.5/k3));

%让0.5以前的都为0
h_vled1(1:length(vled(vled<0))) = 0;
h_vled2(1:length(vled(vled<0))) = 0;
h_vled3(1:length(vled(vled<0))) = 0;

%%
%画图函数
figure(1);
grid on;
hold on;

%k = 2
plot(vled,h_vled1,'R','Linewidth',0.8);
hold on;

%k = 3
plot(vled,h_vled2,'G','Linewidth',0.8);
hold on;

%k = 50
plot(vled,h_vled3,'B','Linewidth',0.8);
hold on;

title('LED-Curve fitting');
xlabel('LED电压');
ylabel('LED电流');
legend('k=2','k=3','k=50');
set(gca,'XMinorGrid','on');
%%
%***********************************第二问****************************%
%%

%拟合
v = [2.76, 2.77, 2.8, 2.85, 2.9, 2.95, 3, 3.05, 3.08, 3.1, 3.2, 3.3, 3.4];
v1 = v-2.76; 
a = [0, 0.00001, 0.015, 0.03, 0.055, 0.09, 0.12, 0.15, 0.19, 0.21, 0.35, 0.475, 0.65];
p = polyfit(v1,a,2);
a_fit = p(1)*v1.^2+p(2)*v1+p(3);

figure(2);
grid on;
hold on;

plot(v,a,'R','Linewidth',0.8);
hold on;

plot(v,a_fit,'G','Linewidth',0.8);
hold on;

title('polyfit fitting');
xlabel('电压');
ylabel('电流');
legend('原始','拟合');
set(gca,'XMinorGrid','on');

%%
%遍历法求参数
% 定义i和k和v的值范围与步长
i_range = 0.01:0.01:3;
k_range = 0.01:0.01:5;
V_range = 0:0.01:5;
% 初始化变量以存储最优值和最小误差
opt_i = 0;
opt_k = 0;
min_error = Inf;  %无穷大

% 嵌套循环以迭代所有i和k的组合
for i = i_range
    for k = k_range
        % 计算当前i和k时拟合值和datasheet均方误差
       f_fit_0 = p(1)*v1.^2+p(2)*v1+p(3);  % 求解二次拟合的代替公式（4）的函数
        f_fit = f_fit_0./((1+(f_fit_0/i).^(2*k)).^(0.5/k));%拟合的最终电流
        error_single = (f_fit-a).^2;
        error = sum(error_single);
        
        % 如果找到新的最小值，则更新最优值和最小误差
        if error < min_error
            min_error = error;
            opt_i = i;
            opt_k = k;
        end
    end
end

% 显示最优值和最小误差
fprintf('最优值：i = %.2f，k = %.2f\n最小误差：%.4f\n', opt_i, opt_k, min_error);
%opt_k = 3;

V_final2 = p(1)*V_range.^2+p(2)*V_range+p(3);
I_final = V_final2./((1+(V_final2/opt_i).^(2*opt_k)).^(0.5/opt_k));

%画图
figure(3);
grid on;
plot(V_range+2.76, I_final,'R','LineWidth',0.8);
title('最优拟合曲线');
xlabel('电压');
ylabel('电流');
legend('拟合曲线 k=2.44 imax=1.7A');
set(gca,'XMinorGrid','on');
set(gca,'YMinorGrid','on');

% clear;