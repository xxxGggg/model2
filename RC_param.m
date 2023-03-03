clear;
clc;
format long;


bp = bodeoptions;
bp.FreqScale = 'linear';
bp.Grid = 'on';
bp.PhaseWrapping = 'on';
bp.XLim=[0.1,130e3 / 2];
bp. FreqUnits = 'Hz';

fsw = 130e3;%开关频率
Ts_Control = 1 / fsw;
V_o = 390;%输出电压参考
L = 125e-6*2;%715.3e-6;%电感大小

z = tf('z',Ts_Control);
s = tf('s');

Td = 1.5 * Ts_Control;

wsb = 2 * pi * fsw;
fsp = 0:1:fsw / 2;
ws = 2 * pi * fsp;

%%电流环设计
phi_m = 60 * pi / 180;%期望裕度
phi_0 = 0;

wc = (pi/2 - phi_0 - phi_m) / Td;%计算截止频率

%phi_1 = 2 * pi * 50 * Td;

Kp = wc * L / V_o;
%Ki = Kp * wc * phi_0 / cos(phi_1);
Gc = Kp;



G_d = exp(-Td * s);
H_op_s = G_d * V_o / (L * s) * Gc;
H_cp_s = feedback(H_op_s,1);

Krmax=[];   a0=0.5;a1=0.25;
    MagTc=wc./sqrt(wc^2+ws.^2-2*wc*ws.*sin(ws*Td)); %%N(w)
    PhaTc=-pi+atan2(ws.*cos(ws*Td),(ws.*sin(ws*Td)-wc));  %%0(w)
%         [MagTc,PhaTc]=bode(Tc,wsu);    
    for m=0:10
            ThetaT=PhaTc+m*ws*Ts_Control;  
            CosThetaT=cos(ThetaT);  %%a(w)
            %bw = 1 ./ (a0 + 2 * a1 * cos(ws * Ts_Control)) ./ (a0 + 2 * a1 * cos(ws * Ts_Control)) - sin(ThetaT) .* sin(ThetaT); %%b(w)
            CoswTs=cos(ws*Ts_Control);  
            fx=1./MagTc.*(CosThetaT+sqrt(CosThetaT.^2+((1/1.0542)^2*4./(1+CoswTs).^2))-1); %n=1,0.9999;n=3.1.0542
            %fx = (CosThetaT + sqrt(bw)) / MagTc; 
            Krmax(m+1)=min(fx);
    end
        
    figure(2);
    m=0:10;
    plot(m,Krmax,'o-','linewidth',1.2);
    xlabel('\it m');
    ylabel('\it K_{rmax}');  
    grid on;    
    axis([0 10 -1 2]);

    Kr =0.7;
    m = 2;
    fg = 50;
    N = fsw / fg;
    Ts = Ts_Control;
    Qs = (0.25+0.5*exp(-Ts*s)+0.25*exp(-2*Ts*s));
    %Gf = exp(m*Ts*s);
    Gr =Kr * Qs * exp((-N+m)*Ts*s) / (1 - Qs * exp((-N+m)*Ts*s)*exp(-m*Ts*s));
    Gic = (1 + Gr) * Kp;

    Giv = 1 / s / L;
    G = V_o / L / s;
    id_ref = 2000 / 220 * sqrt(2);
    Kf = 1 / 220 / sqrt(2);
    Y1 = Giv / (1 + Gic * G * G_d);
    Y2 = Kf * id_ref * (Gic * G_d * G) / (1 + Gic * G * G_d);
    Y = Y1 + Y2;
   
    [reY,imY] = nyquist(Y,ws);

    figure(4)
    plot(ws / 2 / pi,reY(1,:),'b');
    xlabel('\itfreq');
    ylabel('Re\{\itY(jω_s)\}');
    Lg = 0.2e-3;
    Yg = 1 / Lg / s;
    figure(5);
    bode(Y,'b',Yg,'--r',bp);










