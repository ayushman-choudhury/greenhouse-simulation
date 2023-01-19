% Ayushman Choudhury - Greenhouse Gases x Daisy World
% This is the main program, 
% which lays out the Daisy World functionality.
% My helper function for Greenhouse Gases is in a separate file

clear variables
modelNum = 8; % <---------- Change this as desired ----------

% Model 0: Base Daisy World
% Model 1: Constant emission
% Model 2: Emission = Absorption
% Model 3: Better scheme, R = 0
% Model 4: Learn slowly, R += 10%
% Model 5: Learn slowly, R += 100%
% Model 6: Learn faster, R += 2%
% Model 7: Learn faster, R += 10%
% Model 8: Learn faster, R += 100%

if modelNum < 2
    initE = 0;
else
    initE = 10;
end

% Pre-defined constants

initCoverGround = 0.9;
initCoverWhite = 0.05;
initCoverBlack = 0.05;

Ag = 0.5;
Aw = 0.8;
Ab = 0.2;

L = 1;
S0 = 1380;
K = 0.6;
sigma = 5.6 * 10^-8;

b = 3.265 * 10^-3;
T0 = 293.7;
D = 0.3;

dt = 1;
R = 0;

% Set up matrix

% I use a matrix to store all my data:
% 1 - Year
% 2 - White Cover (marching equations)
% 3 - Black Cover (marching equations)
% 4 - Barren Ground Cover (subtract)
% 5 - Overall Albedo (cover * albedo per surface type)
% 6 - Atmospheric temperature (depends on Overall Albedo)
% 7 - Surface temperature (derived from Atmospheric)
% 8 - White temperature (Surface Temp & Albedo)
% 9 - Black temperature (same)
% 10 - White birth rate (depends on white temp)
% 11 - Black birth rate (depends on black temp)
% 12 - Greenhouse Effect 

% Load in initial conditions from given constants
years = 100;
data = zeros(12, years + 1);
data(2,1) = initCoverWhite;
data(3,1) = initCoverBlack;
data(12,1) = initE;

% Step for 100 years

for iter = 1:(years + 1)
    year = iter - 1;
    data(1,iter) = year;

    % Computations
    Cw = data(2,iter);
    Cb = data(3,iter);
    Cg = 1 - Cw - Cb;
    A = (Cg * Ag) + (Cw * Aw) + (Cb * Ab);

    % New: greenhouse gases
    E_old = data(12,iter);
    if iter < 2
        [E, R] = greenhouse_effect_final(modelNum, E_old, dt, iter, ...
        Cw, Cw, Cb, Aw, Ab, R);
    else
        [E, R] = greenhouse_effect_final(modelNum, E_old, dt, iter, ...
        data(2,iter - 1), Cw, Cb, Aw, Ab, R);
    end


    % New Ta calculation - all else is same
    Ta = ((L * S0 * (1 - A)) / (4 * sigma))^0.25 + E;
    Ts = (2 * Ta^4)^0.25;
    Tw = ((1 - K) * ((L * S0) / (4 * sigma)) * (A - Aw) + Ts^4)^0.25;
    Tb = ((1 - K) * ((L * S0) / (4 * sigma)) * (A - Ab) + Ts^4)^0.25;
    bw = 1 - b * (T0 - Tw)^2;
    bb = 1 - b * (T0 - Tb)^2;

    % Make sure birth rate doesn't become negative
    if bb < 0
        bb = 0;
    end
    if bw < 0
        bw = 0;
    end

    % Store current values in matrix
    data(4,iter) = Cg;
    data(5,iter) = A;
    data(6,iter) = Ta;
    data(7,iter) = Ts;
    data(8,iter) = Tw;
    data(9,iter) = Tb;
    data(10,iter) = bw;
    data(11,iter) = bb;

    if iter < (years + 1)
        % If we're not in the last column,
        % put the next Cw and Cb values in the matrix
        
        % For the homework, I had used a Forward Euler
        %next_Cw = Cw * (1 + dt*bw*(1 - Cw - Cb) - dt * D);
        %next_Cb = Cb * (1 + dt*bb*(1 - Cb - Cw) - dt * D);

        % New: Runge-Kutta 2
        wK1 = dt * (bw * Cw * (1 - Cw - Cb) - D * Cw);
        Cw_mid = Cw + wK1;
        wK2 = dt * (bw * Cw_mid * (1 - Cw_mid - Cb) - D * Cw_mid);
        Cw_next = Cw + (wK1 + wK2) / 2;

        bK1 = dt * (bb * Cb * (1 - Cb - Cw) - D * Cb);
        Cb_mid = Cb + bK1;
        bK2 = dt * (bb * Cb_mid * (1 - Cb_mid - Cw) - D * Cb_mid);
        Cb_next = Cb + (bK1 + bK2) / 2;

        data(2,iter + 1) = Cw_next;
        data(3,iter + 1) = Cb_next;
        data(12,iter + 1) = E;
    end
end

%Plots

yrs = data(1,:);

subplot(2,1,1);
plot(yrs, data(6,:));
title('Atmospheric Temperature over time');
xlabel('Time in Years');
ylabel('Atmos. Temp.');
legend({'Atmospheric Temperature'}, 'Location', 'northeast');

subplot(2,1,2);
plot(yrs, data(2,:), 'c', yrs, data(3,:), 'k', yrs, data(4,:), 'r');
title('Cover Fractions by surface types over time');
xlabel('Time in Years');
ylabel('Cover Fractions');
legend({'White Daisies', 'Black Daisies', 'Barren Ground'}, ...
    'Location', 'northeast');
