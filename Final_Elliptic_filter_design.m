% clc;
% clear;
% close all;
M = input("Enter the Filter Number : ");
delta1 = input(" Enter the value of Passband Tolerance : ");
delta2 = input(" Enter the value of Stopband Tolerance : ");
Tw = input(" Enter the value of transition width : ");
fs = input("Enter the value of Sampling Frequency : ")
Q = floor(M/11)  % Quotient 
R = rem(M,11)    % Remainder
%% Filter 1 Calculations
% Calculation of Stopband Range for Filter 1
fs1 = 40 + 5*Q
fs2 = 70 + 5*Q
fprintf('Stopband Range of filter 1 should be: %.2f kHz to %.2f kHz\n', fs1, fs2);
% Calculation of Passband Range for Filter 1
fp1 = fs1 - Tw
fp2 = fs2 + Tw
fprintf('Stopband Range of filter 1 should be: 0 to %.2f kHz and %.2f kHz to inf\n', fp1, fp2);

%% Bilinear Transformation - tan(pi*f/fs) for Filter 1
omegas1 = tan(pi*fs1/fs)
omegas2 = tan(pi*fs2/fs)
omegap1 = tan(pi*fp1/fs)
omegap2 = tan(pi*fp2/fs)
B = omegap2 - omegap1
omega0 = sqrt(omegap2*omegap1)
fprintf('Stopband Range of filter 1 should be: %.2f kHz to %.2f kHz\n', omegas1, omegas2);
fprintf('Stopband Range of filter 1 should be: 0 to %.2f kHz and %.2f kHz to inf\n', omegap1, omegap2);

%% Mapping of BSP to Butterworth LPF - B*Omega/(omega0^2 - omega^2) for Filter 1
omegaLp1 = B*omegap1/(omega0^2 - omegap1^2)
omegaLp2 = B*omegap2/(omega0^2 - omegap2^2)
omegaLs1 = B*omegas1/(omega0^2 - omegas1^2)
omegaLs2 = B*omegas2/(omega0^2 - omegas2^2)
fprintf('Stopband Range of filter 1 should be: %.2f kHz to %.2f kHz\n', omegaLs1, omegaLs2);
fprintf('Stopband Range of filter 1 should be: 0 to %.2f kHz and %.2f kHz to inf\n', omegaLp1, omegaLp2);
omegaLs = min(abs(omegaLs1),abs(omegaLs2))
% Calculation of order
k = 1/omegaLs
k_ = sqrt(1 - k^2)
epsilon = sqrt(2*delta1 - delta1^2)/sqrt(1 - 2*delta1 + delta1^2)
k1 = epsilon/sqrt(1/delta2^2 -1)
k1_ = sqrt(1 - k1^2)
fprintf("Value of k : ")
disp(k)
fprintf("Value of kdass : ")
disp(k_)
fprintf("Value of k1 : ")
disp(k1)
fprintf("Value of k1dass : ")
disp(k1_)
%integral calculation
f1 = @(y) 1 ./ sqrt(1 - k^2 * sin(y).^2);
K = integral(f1,0,pi/2)
disp(K)
f2 = @(y) 1 ./ sqrt(1 - k_^2 * sin(y).^2);
K_ = integral(f2,0,pi/2)
disp(K_)
f3 = @(y) 1 ./ sqrt(1 - k1^2 * sin(y).^2);
K1 = integral(f3,0,pi/2)
disp(K1)
f4 = @(y) 1 ./ sqrt(1 - k1_^2 * sin(y).^2);
K1_ = integral(f4,0,pi/2)
disp(K1_)
f5 = inv_sc_elliptic(1/epsilon,k1_)
N = (K*K1_)/(K_*K1)
N = ceil(N)
fprintf(" Order of Filter is : ")
disp(N)

zeros = [];  % Initialize

if mod(N,2) == 1
    % --- For odd N ---
    % o = 0, 2, 4, ..., (N-1)
    o_vals = 0:2:(N-1);
else
    % --- For even N ---
    % o = 1, 3, 5, ..., (N-1)
    o_vals = 1:2:(N-1);
end

for o = o_vals
    [sn, cn, dn] = ellipj_from_scratch(o * K / N, k^2);
    z_val = 1j ./ (k .* sn);       % Zero locations
    zeros = [zeros, z_val, -z_val]; % Add conjugate pair
end

fprintf('\nZeros:\n');
disp(zeros.');

%poles
f5 = inv_sc_elliptic(1/epsilon,k1_^2)
v = (K*f5)/(N*K1)

poles = [];

if mod(N,2) == 1
    % N is odd
    l_vals = 0:2:(N-1);
else
    % N is even
    l_vals = 1:2:(N-1);
end

for l = l_vals
    [sn, cn, dn] = ellipj_from_scratch(l * K / N, k^2);
    [sn_, cn_, dn_] = ellipj_from_scratch(v, k_^2);

    % Pole formula
    p_val = ((-1)*cn * dn * sn_ * cn_ + 1j * sn * dn_) / (1 - dn^2 * sn_^2);

    % Add only one copy if real, else add conjugate pair
    if abs(imag(p_val)) < 1e-12
        poles = [poles, p_val];   % real pole → add once
    else
        poles = [poles, p_val, conj(p_val)];  % complex pole pair
    end
end

fprintf('\nPoles:\n');
disp(poles.');
%% -------------------- Transfer Function H(s) --------------------
[num_lp, den_lp] = zp2tf(zeros.', poles.', 1);
Hlp = tf(num_lp, den_lp);
H0q = abs(polyval(num_lp,0) / polyval(den_lp,0));  % DC gain
num_scaled = num_lp / H0q;                 % scale numerator for |H(0)| = 1
Hlp = tf(num_scaled, den_lp);
disp('Butterworth Lowpass Transfer Function H(s) of filter 1:');
Hlp

%% ----Low pass Elliptic filter -----
figure;
bode(Hlp);
grid on;
title('Magnitude and Phase Response of Elliptic Filter H(s)');

%% ----- Pole-Zero Plot -----
figure;
pzmap(Hlp);
grid on;
title('Pole-Zero Map of Elliptic Filter H(s)');

%% ----- Frequency Response -----
figure;
freqs(num_scaled, den_lp);
title('Frequency Response |H(j\omega)| of Elliptic Filter');
xlabel('Frequency (rad/s)');
ylabel('Magnitude');
grid on;

%% ---------- Bandstop Substitution s -> (B*s)/(s^2+Ω0^2) ----------
syms s z
qsubs_expr = (B*s)/(s^2 + omega0^2); % for First Filter
num_bp = poly2sym(num_scaled,s);
den_bp = poly2sym(den_lp,s);

qnum_bp_sub = expand(subs(num_bp, s, qsubs_expr));
qden_bp_sub = expand(subs(den_bp, s, qsubs_expr));
qHbp = qnum_bp_sub/qden_bp_sub;
fprintf('Analog Bandpass Transfer Function (First Filter) qHbp(s)\n');
qHbp
%% ---------- Extract Analog Bandstop Coefficients ----------
% Convert symbolic transfer function to numerator and denominator polynomials
[num_bs_analog, den_bs_analog] = numden(qHbp);

% Extract coefficients
num_bs_analog_coeff = sym2poly(num_bs_analog);
den_bs_analog_coeff = sym2poly(den_bs_analog);

% Normalize by leading denominator coefficient
k_analog = den_bs_analog_coeff(1);
num_bs_analog_coeff = num_bs_analog_coeff / k_analog;
den_bs_analog_coeff = den_bs_analog_coeff / k_analog;

fprintf('\n=== FILTER 1 ANALOG BANDSTOP COEFFICIENTS ===\n');
fprintf('Analog Numerator coefficients (s-domain): ');
fprintf('%.6f ', num_bs_analog_coeff);
fprintf('\n');
fprintf('Analog Denominator coefficients (s-domain): ');
fprintf('%.6f ', den_bs_analog_coeff);
fprintf('\n');

%% ---------- Bilinear Transformation (s -> (z-1)/(z+1)) ---------
qHbp_z = simplify(subs(qHbp, s, (z-1)/(z+1)));

disp('Discrete Transfer Function (First Filter) H(z):');
pretty(qHbp_z);

%% ---------- Convert to numeric coefficients ----------
[nz1, dz1] = numden(qHbp_z);
nz1 = sym2poly(expand(nz1));
dz1 = sym2poly(expand(dz1));
% --- Normalize by leading denominator coefficient (you already do this)
k12  = dz1(1);
nz1 = nz1 / k12;
dz1 = dz1 / k12;
disp('Numerator coefficients of analog BSF-I in z:');
disp(mat2str(nz1));
disp('Denominator coefficients of analog BSF-I in z:');
disp(mat2str(dz1));
fvtool(nz1,dz1)

%% Filter 2 Calculations
% Calculation of Stopband Range for Filter 2
fs11 = 170 + 5*R
fs22 = 200 + 5*R
fprintf('Stopband Range of filter 2 should be: %.2f kHz to %.2f kHz\n', fs11, fs22);
% Calculation of Passband Range for Filter 2
fp11 = fs11 - Tw
fp22 = fs22 + Tw
fprintf('Passband Range of filter 2 should be: 0 to %.2f kHz and %.2f kHz to inf\n', fp11, fp22);

%% Bilinear Transformation - tan(pi*f/fs) for Filter 2
omegas11 = tan(pi*fs11/fs)
omegas22 = tan(pi*fs22/fs)
omegap11 = tan(pi*fp11/fs)
omegap22 = tan(pi*fp22/fs)
B1 = omegap22 - omegap11
omega01 = sqrt(omegap22*omegap11)
fprintf('Stopband Range of filter 2 should be: %.2f kHz to %.2f kHz\n', omegas11, omegas22);
fprintf('Stopband Range of filter 2 should be: 0 to %.2f kHz and %.2f kHz to inf\n', omegap11, omegap22);

%% Mapping of BSP to Butterworth LPF - B*Omega/(omega0^2 - omega^2) for Filter 2
omegaLp11 = B1*omegap11/(omega01^2 - omegap11^2)
omegaLp22 = B1*omegap22/(omega01^2 - omegap22^2)
omegaLs11 = B1*omegas11/(omega01^2 - omegas11^2)
omegaLs22 = B1*omegas22/(omega01^2 - omegas22^2)
fprintf('Stopband Range of filter 2 should be: %.2f kHz to %.2f kHz\n', omegaLs11, omegaLs22);
fprintf('Stopband Range of filter 2 should be: 0 to %.2f kHz and %.2f kHz to inf\n', omegaLp11, omegaLp22);
omegaLs1 = min(abs(omegaLs11),abs(omegaLs22))
% Calculation of order
k = 1/omegaLs1
k_ = sqrt(1 - k^2)
epsilon = sqrt(2*delta1 - delta1^2)/sqrt(1 - 2*delta1 + delta1^2)
k1 = epsilon/sqrt(1/delta2^2 -1)
k1_ = sqrt(1 - k1^2)
fprintf("Value of k : ")
disp(k)
fprintf("Value of kdass : ")
disp(k_)
fprintf("Value of k1 : ")
disp(k1)
fprintf("Value of k1dass : ")
disp(k1_)
%integral calculation
f1 = @(y) 1 ./ sqrt(1 - k^2 * sin(y).^2);
K = integral(f1,0,pi/2)
disp(K)
f2 = @(y) 1 ./ sqrt(1 - k_^2 * sin(y).^2);
K_ = integral(f2,0,pi/2)
disp(K_)
f3 = @(y) 1 ./ sqrt(1 - k1^2 * sin(y).^2);
K1 = integral(f3,0,pi/2)
disp(K1)
f4 = @(y) 1 ./ sqrt(1 - k1_^2 * sin(y).^2);
K1_ = integral(f4,0,pi/2)
disp(K1_)
f5 = inv_sc_elliptic(1/epsilon,k1_)
N2 = (K*K1_)/(K_*K1)
N2 = ceil(N2)
fprintf(" Order of 2nd Filter is : ")
disp(N2)
% Calculation of Poles
zeros1 = [];  % Initialize

if mod(N,2) == 1
    % --- For odd N ---
    % o = 0, 2, 4, ..., (N-1)
    o_vals1 = 0:2:(N2-1);
else
    % --- For even N ---
    % o = 1, 3, 5, ..., (N-1)
    o_vals1 = 1:2:(N2-1);
end

for o = o_vals1
    [sn1, cn1, dn1] = ellipj_from_scratch(o * K / N2, k^2);
    z_val1 = 1j ./ (k .* sn1);       % Zero locations
    zeros1 = [zeros1, z_val1, -z_val1]; % Add conjugate pair
end

fprintf('\nZeros:\n');
disp(zeros1.');

%Calculation of poles
f5 = inv_sc_elliptic(1/epsilon,k1_^2)
v1 = (K*f5)/(N2*K1)

poles1 = [];

if mod(N2,2) == 1
    % N is odd
    l_vals = 0:2:(N2-1);
else
    % N is even
    l_vals = 1:2:(N2-1);
end

for l = l_vals
    [sn1, cn1, dn1] = ellipj_from_scratch(l * K / N2, k^2);
    [sn_, cn_, dn_] = ellipj_from_scratch(v, k_^2);

    % Pole formula
    p_val = ((-1)*cn1 * dn1 * sn_ * cn_ + 1j * sn1 * dn_) / (1 - dn1^2 * sn_^2);

    % Add only one copy if real, else add conjugate pair
    if abs(imag(p_val)) < 1e-12
        poles1 = [poles1, p_val];   
    else
        poles1 = [poles1, p_val, conj(p_val)];  % complex pole pair
    end
end

fprintf('\nPoles:\n');
disp(poles1.');
%% -------------------- Transfer Function H(s) --------------------
[num_lp, den_lp] = zp2tf(zeros1.', poles1.', 1);
Hlp = tf(num_lp, den_lp);
H0q = abs(polyval(num_lp,0) / polyval(den_lp,0));  % DC gain
num_scaled = num_lp / H0q;                 % scale numerator for |H(0)| = 1
Hlp = tf(num_scaled, den_lp);
disp('Butterworth Lowpass Transfer Function H(s) of filter 2:');
Hlp

%% ----Low pass Elliptic filter -----
figure;
bode(Hlp);
grid on;
title('Magnitude and Phase Response of Elliptic Filter H(s)');

%% ----- Pole-Zero Plot -----
figure;
pzmap(Hlp);
grid on;
title('Pole-Zero Map of Elliptic Filter H(s)');

%% ----- Frequency Response -----
figure;
freqs(num_scaled, den_lp);
title('Frequency Response |H(j\omega)| of Elliptic Filter');
xlabel('Frequency (rad/s)');
ylabel('Magnitude');
grid on;

%% ---------- Bandstop Substitution s -> (B*s)/(s^2+Ω0^2) ----------
syms s z
qsubs_expr = (B1*s)/(s^2 + omega01^2); 
num_bp = poly2sym(num_scaled,s);
den_bp = poly2sym(den_lp,s);

qnum_bp_sub = expand(subs(num_bp, s, qsubs_expr));
qden_bp_sub = expand(subs(den_bp, s, qsubs_expr));
qHbp = qnum_bp_sub/qden_bp_sub;
fprintf('Analog Bandpass Transfer Function (2nd Filter) qHbp(s)\n');
qHbp;
%% ---------- Extract Analog Bandstop Coefficients ----------
% Convert symbolic transfer function to numerator and denominator polynomials
[num_bs_analog, den_bs_analog] = numden(qHbp);

% Extract coefficients
num_bs_analog_coeff = sym2poly(num_bs_analog);
den_bs_analog_coeff = sym2poly(den_bs_analog);

% Normalize by leading denominator coefficient
k_analog = den_bs_analog_coeff(1);
num_bs_analog_coeff = num_bs_analog_coeff / k_analog;
den_bs_analog_coeff = den_bs_analog_coeff / k_analog;

fprintf('\n=== FILTER 1 ANALOG BANDSTOP COEFFICIENTS ===\n');
fprintf('Analog Numerator coefficients (s-domain): ');
fprintf('%.6f ', num_bs_analog_coeff);
fprintf('\n');
fprintf('Analog Denominator coefficients (s-domain): ');
fprintf('%.6f ', den_bs_analog_coeff);
fprintf('\n');
%% ---------- Bilinear Transformation (s -> (z-1)/(z+1)) ---------
qHbp_z = simplify(subs(qHbp, s, (z-1)/(z+1)));

disp('Discrete Transfer Function (2nd Filter) H(z):');
pretty(qHbp_z);

%% ---------- Convert to numeric coefficients ----------
[nz11, dz11] = numden(qHbp_z);
nz11 = sym2poly(expand(nz11));
dz11 = sym2poly(expand(dz11));
% --- Normalize by leading denominator coefficient
k11  = dz11(1);
nz11 = nz11 / k11;
dz11 = dz11 / k11;
disp('Numerator coefficients of analog BSF-II in z:');
disp(mat2str(nz11));
disp('Denominator coefficients of analog BSF-II in z:');
disp(mat2str(dz11));
fvtool(nz11,dz11)

%Cascade of 2 Band stop Filter
x = conv(nz1,nz11)
y = conv(dz1,dz11)
disp('Numerator coefficients of analog Multi - BSF in z:');
disp(mat2str(x));
disp('Denominator coefficients of analog Multi - BSF in z:');
disp(mat2str(y));
fvtool(x,y)