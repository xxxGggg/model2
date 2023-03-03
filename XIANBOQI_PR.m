
clear;
clc;
%bode图选项
bp = bodeoptions;
bp.XLim=[0.1,130e3];
bp. FreqUnits = 'Hz';

fsw = 130e3;%开关频率
Vo = 390;%输出电压参考
L = 125e-6*2;%715.3e-6;%电感大小

s = tf('s');

%%电流环设计
phi_m = 60 * pi / 180;%期望裕度
phi_0 = 0.1;
Ts_Control = 1 / fsw ; %控制周期
Td = 1.5 * Ts_Control;%延时
w_c = (pi/2 - phi_0 - phi_m) / Td;%计算截止频率

phi_1 = 2 * pi * 50 * Td * 0;

Kp = w_c * L / Vo;
Ki = Kp * w_c * phi_0 / cos(phi_1);

omega1 = 2 * pi * 100;
wrc = 31.4;
Gic = 1 - ( (s*cos(phi_1) - omega1 * sin(phi_1)) / (s^2 + wrc * s + omega1^2));




figure(1);
bode(Gic,bp);



