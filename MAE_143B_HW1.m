% MAE 143B HW 1

%% Problem 1

s = tf('s');

wg = 10;
alpha = 14; % estimated from Figure 10.17c

z = wg/sqrt(alpha);
p = wg*sqrt(alpha);

Dlead = (s + z)/(s + p);

figure
bode(Dlead)
grid on
title('Bode Plot of Lead Compensator D_{lead}(s)')

%% Problem 2a


s = tf('s');

z2a = 0.885;
p2a = 0.00885;

Dlag = (s + z2a)/(s + p2a);

figure
bode(Dlag)
grid on
title('Bode Plot of Lag Compensator D_{lag}(s)')

%% Problem 2b

s = tf('s');

z2b = 0.485;
p2b = 0.0485;

Ddoublelag = ((s + z2b)/(s + p2b))^2;

figure
bode(Dlag, Ddoublelag)
grid on
legend('D_{lag}(s)','D_{double-lag}(s)')
title('Bode Plot of Lag and Double-Lag Compensators')

%% Problem 3a

wg = 10;
phase_desired = -5;

wc = 100;       % initial guess

for k = 1:10000 % iteration to find wc

    DLPF4 = RR_LPF_butterworth(4,wc);

    Dwg = RR_evaluate(DLPF4,1i*wg);

    phase = angle(Dwg)*180/pi;

    if abs(phase - phase_desired) < 0.001
        break
    end

    if phase < phase_desired
        wc = wc + 0.1;
    else
        wc = wc - 0.1;
    end
end

figure
RR_bode(DLPF4)
title('4th-Order Butterworth Low-Pass Filter')

%% Problem 3b

wg = 10;
delta = 0.001;
phase_desired = -5;

wc = 500;       % initial guess

for k = 1:10000

    Dinv = RR_LPF_inv_chebyshev(4,delta,wc);

    Dwg = RR_evaluate(Dinv,1i*wg);

    phase = angle(Dwg)*180/pi;

    if abs(phase - phase_desired) < 0.001
        break
    end

    if phase < phase_desired
        wc = wc + 0.1;
    else
        wc = wc - 0.1;
    end
end

figure
RR_bode(DLPF4)
RR_bode(Dinv)
legend('Butterworth','Inverse Chebyshev')
title('4th-Order Inverse Chebyshev Filter')

Rs = -20*log10(delta);

[num,den] = cheby2(4,Rs,wc,'s');

Dinv = tf(num,den)

%% Problem 4b

s = tf('s');

z_lead = 2.672;
p_lead = 37.416;

D_lead = (s + z_lead)/(s + p_lead);

z_lag = 0.485;
p_lag = 0.0485;

D_doublelag = ((s + z_lag)/(s + p_lag))^2;

num_inv = [0.001 0 7666 0 7.346e9];
den_inv = [1 758.1 2.874e5 6.41e7 7.346e9];

D_inv = tf(num_inv,den_inv);

% CT loop-shaping controller without K
Dloop_noK = D_lead * D_doublelag * D_inv;

disp('Continuous-time controller without K:')
Dloop_noK


% Tustin with prewarping

h = 0.001;
wg = 10;

opts = c2dOptions('Method','tustin', ...
    'PrewarpFrequency',wg);

Dz = c2d(Dloop_noK,h,opts);

disp('Discrete-time controller without K:')
Dz


% Extract numerator and denominator coefficients

[numz,denz] = tfdata(Dz,'v');

% Normalize so first denominator coefficient = 1
numz = numz/denz(1);
denz = denz/denz(1);

% Difference equation coefficients
b = numz;
a = denz(2:end);


% Print to 4 significant digits

fprintf('\n--- a coefficients ---\n')
for j = 1:length(a)
    fprintf('a%d = %.4g\n',j,a(j));
end

fprintf('\n--- b coefficients ---\n')
for j = 1:length(b)
    fprintf('b%d = %.4g\n',j-1,b(j));
end

%% Problem 5a

G = RR_tf(100,[1 0 -100]);

z = 10;
p = 20;

Dsimple0 = RR_tf([1 z],[1 p]);
Dsimple0

g.K = logspace(-2,2,1000);
g.axes = [-30 15 -20 20];

figure
RR_rlocus(G,Dsimple0,g)

hold on

% Desired poles from tr = 0.18 and Mp = 15%
sd1 = -5 + 1i*8.6603;
sd2 = -5 - 1i*8.6603;

plot(real(sd1),imag(sd1),'rx','MarkerSize',10,'LineWidth',2)
plot(real(sd2),imag(sd2),'rx','MarkerSize',10,'LineWidth',2)

title('Root Locus for D_{simple}(s)')

%% Problem 5b
s = tf('s');

G = 100/(s^2 - 100);

Dsimple = 3*(s + 10)/(s + 20);

Dlead = (s + 2.672)/(s + 37.416);

Ddoublelag = ((s + 0.485)/(s + 0.0485))^2;

num_inv = [0.001 0 7666 0 7.346e9];
den_inv = [1 758.1 2.874e5 6.41e7 7.346e9];

Dinv = tf(num_inv,den_inv);

Dloop_noK = Dlead*Ddoublelag*Dinv;

wg = 10;

L0_wg = evalfr(G*Dloop_noK,1i*wg);

K = 1/abs(L0_wg);

fprintf('K = %.6f\n',K)

Dloop = K*Dloop_noK;

Lsimple = G*Dsimple;
Lloop = G*Dloop;

figure
rlocus(G*((s+10)/(s+20)))
grid on
title('Root Locus - Simple Lead Controller')

figure
rlocus(G*Dloop_noK)
grid on
title('Root Locus - Loop Shaping Controller')

figure
bode(Lsimple,Lloop)
grid on
legend('Simple Lead','Loop Shaping')
title('Open-Loop Bode Plots')

Tsimple = feedback(Lsimple,1);
Tloop = feedback(Lloop,1);

figure
step(Tsimple,Tloop)
grid on
legend('Simple Lead','Loop Shaping')
title('Closed-Loop Step Responses')

[mag,phase] = bode(Lloop,wg);

mag = squeeze(mag);
phase = squeeze(phase);

fprintf('Loop-shaping magnitude at wg = %.6f\n',mag)
fprintf('Loop-shaping phase at wg = %.4f deg\n',phase)
fprintf('Approximate phase margin = %.4f deg\n',180 + phase)

info_simple = stepinfo(Tsimple);
info_loop = stepinfo(Tloop);

disp('Simple Lead Step Response:')
disp(info_simple)

disp('Loop Shaping Step Response:')
disp(info_loop)