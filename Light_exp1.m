
%%
%�����Ļ�����ʽ�����
R = 1;
imax = 0.5;
k1 = 2;k2 = 3;k3 = 50;
vled = -0.5:0.1:1;

%����Ҫ��.^,������ʱ��Ҳ�ǣ��������г���ʱҪ�õ�
%��ʽ�е�h
f_vled = vled/R;
h_vled1 = f_vled./((1+(f_vled/imax).^(2*k1)).^(0.5/k1));
h_vled2 = f_vled./((1+(f_vled/imax).^(2*k2)).^(0.5/k2));
h_vled3 = f_vled./((1+(f_vled/imax).^(2*k3)).^(0.5/k3));

%��0.5��ǰ�Ķ�Ϊ0
h_vled1(1:length(vled(vled<0))) = 0;
h_vled2(1:length(vled(vled<0))) = 0;
h_vled3(1:length(vled(vled<0))) = 0;

%%
%��ͼ����
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
xlabel('LED��ѹ');
ylabel('LED����');
legend('k=2','k=3','k=50');
set(gca,'XMinorGrid','on');
%%
%***********************************�ڶ���****************************%
%%

%���
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
xlabel('��ѹ');
ylabel('����');
legend('ԭʼ','���');
set(gca,'XMinorGrid','on');

%%
%�����������
% ����i��k��v��ֵ��Χ�벽��
i_range = 0.01:0.01:3;
k_range = 0.01:0.01:5;
V_range = 0:0.01:5;
% ��ʼ�������Դ洢����ֵ����С���
opt_i = 0;
opt_k = 0;
min_error = Inf;  %�����

% Ƕ��ѭ���Ե�������i��k�����
for i = i_range
    for k = k_range
        % ���㵱ǰi��kʱ���ֵ��datasheet�������
       f_fit_0 = p(1)*v1.^2+p(2)*v1+p(3);  % ��������ϵĴ��湫ʽ��4���ĺ���
        f_fit = f_fit_0./((1+(f_fit_0/i).^(2*k)).^(0.5/k));%��ϵ����յ���
        error_single = (f_fit-a).^2;
        error = sum(error_single);
        
        % ����ҵ��µ���Сֵ�����������ֵ����С���
        if error < min_error
            min_error = error;
            opt_i = i;
            opt_k = k;
        end
    end
end

% ��ʾ����ֵ����С���
fprintf('����ֵ��i = %.2f��k = %.2f\n��С��%.4f\n', opt_i, opt_k, min_error);
%opt_k = 3;

V_final2 = p(1)*V_range.^2+p(2)*V_range+p(3);
I_final = V_final2./((1+(V_final2/opt_i).^(2*opt_k)).^(0.5/opt_k));

%��ͼ
figure(3);
grid on;
plot(V_range+2.76, I_final,'R','LineWidth',0.8);
title('�����������');
xlabel('��ѹ');
ylabel('����');
legend('������� k=2.44 imax=1.7A');
set(gca,'XMinorGrid','on');
set(gca,'YMinorGrid','on');

% clear;