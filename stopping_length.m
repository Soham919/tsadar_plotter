clc
clear
close all

v1 = 10^6;                 % projectile velocity (N)
m1 = 14 * 1.6726e-27;     % Ion mass
m2 = 9.109e-31;           % Electron mass
q1 = 1.602e-19;           % Ion charge
q2 = -1.602e-19;          % Electron charge
n = 5e18;                 % Background plasma density
n2 = n;                    % Same for scatterers
K0 = m1*v1^2/2/(10^6*q1);              % Initial kinetic energy in joules
Te = 500 * 1.602e-19;    % Electron temperature in J 
L = compute_stopping_length(K0, Te, m1, m2, q1, q2, n, n2);
fprintf('Stopping length = %.4e meters\n', L);




function stopping_length = compute_stopping_length(K0, Te, m1, m2, q1, q2, n, n2)
    dK_dl_fun = @(K) stopping_power(K, Te, m1, m2, q1, q2, n, n2); % rate of energy loss
    K = K0;
    dl = 1e-7;          % dx
    stopping_length = 0;

    % Integrate until K ~ 0
    while K > 1e-15  % Avoid numerical noise
        dK = dK_dl_fun(K);
        dK = max(dK, 1e-30);  % Prevent negative or zero stopping power
        dK_actual = dK * dl;  % Energy lost in this step
        K = K - dK_actual;
        stopping_length = stopping_length + dl;

        % Safety break to avoid infinite loop (in case of bugs)
        if stopping_length > 10
            warning('Exceeded 10 meters — check parameters.');
            break;
        end
    end
end


function dK_dl = stopping_power(K, Te, m1, m2, q1, q2, n, n2)
    epsilon0 = 8.854187817e-12;
    e = 1.602176634e-19;
    pi_val = 3.14;
    mr = (m1 * m2) / (m1 + m2);
    v1 = sqrt(2 * K / m1);
    b90 = abs(q1 * q2) / (4 * pi_val * epsilon0 * mr * v1^2);
    lambdaD = sqrt(epsilon0 * Te / (n * e^2));
    Lambda = lambdaD / b90;
    mass_factor = (m1 * m2) / (m1 + m2)^2;
    dK_dl = K * n2 * mass_factor * 8 * pi_val * b90^2 * log(abs(Lambda));
end


