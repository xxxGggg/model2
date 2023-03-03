clc;
clear;




fsu = 0.1:0.1:65000;
ws = 2 * pi * fsu;
s = 1i*ws;

%%EMI阻抗
Lg = 50e-6;
L = 250e-6;
R = 1;
C1 = 2 * 2.2e-6;
C2 = 2 * 0.47e-6;
Z1 = 1 ./ C1 ./ s;
Z2 = L .* s;
Z3 = Lg .* s;
Z4 = 1 ./ C2 ./ s;
Z3 = Z3 .* Z4 ./ (Z3 + Z4);
Z3 = Z3 + Z2;
Z = (Z3 .* Z1) ./ (Z3 + Z1);
Z = 1 ./ Z;
MagZ = abs(Z);
PhaZ = angle(Z) * 180 / pi;
%%PFC阻抗
fsu_PFC = 0.1:0.1:65000;
ws_PFC = 2 * pi * fsu;
s_PFC = 1i*ws;
Voe = 390;
L_PFC = 250e-6;
beta_de_PFC = 220 * sqrt(2) / 390;
R_PFC = 75;
C_PFC = 620e-6;

G1 = -Voe ./ (L_PFC .* s);
G2 = 0;% = -beta_de ./ (L .* s);
G3 = 0;%0.5 * R * beta_de ./ (R .* C .* s + 1);

G = G1 + G2 .* G3;

kcp = -0.0623;%-0.018;%-1.7662e3;
kci = -1.871e3;%-102.2042*5;%-0.0749;
Gic = kcp + kci ./ s;
Giv = 1 ./ s / L_PFC;
id_ref = 2000 / 220 * sqrt(2);
Ts_Control = 1 / 130e3;
Gd = exp(-1.5 * Ts_Control * s);
Kf = 1 / 220 / sqrt(2);
%Y1 = Kf * id_ref * (Gic .* Gd .* G) ./ (1 + Gic .* Gd .* G);
Y1 = Giv ./ (1 + Gic .* G .*Gd );
%Y2 = (R*C.*s + 1) ./  (R*L*C.*s.*s+L*s+R*beta_de*beta_de);
Y2 = 1 ./ s / L_PFC ;
Y3 = Kf * id_ref * (Gic .* Gd .* G) ./ (1 + Gic .* Gd .* G);
Y = Y1 + Y3 - Kf .* Gd .* G ./ (1 + Gic .* Gd .* G);

        
        Z_PFC = Y;
        MagZ_PFC = abs(Z_PFC);
        PhaZ_PFC = angle(Z_PFC) * 180 / pi;



ZN = Z ./ Z_PFC;
Mag_ZN = abs(Z_PFC);
Pha_ZN = angle(ZN);
%figure(1);
%polarplot(Pha_ZN,Mag_ZN);

figure(71);
subplot(2,1,1);
%plot(harmor,jieguo(:,2),'r*');
%hold on;

plot(fsu,20 * log10(MagZ_PFC),'g','LineWidth',2);
hold on;

plot(fsu,20 * log10(MagZ),'b','LineWidth',2);
hold off;

grid on;
subplot(2,1,2);
%plot(harmor,jieguo(:,3),'r*');
%hold on;
plot(fsu,PhaZ_PFC,'g','LineWidth',2);
hold on;
plot(fsu,PhaZ,'b','LineWidth',2);
grid on;
hold off;















