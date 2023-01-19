% Ayushman Choudhury - Greenhouse Gases x Daisy World
% This is the helper function,
% which describes the Greenhouse Effect calculations
% It is separated by Model Number:

% Model 0: Base Daisy World
% Model 1: Constant emission
% Model 2: Emission = Absorption
% Model 3: Better scheme, R = 0
% Model 4: Learn slowly, R += 10%
% Model 5: Learn slowly, R += 100%
% Model 6: Learn faster, R += 2%
% Model 7: Learn faster, R += 10%
% Model 8: Learn faster, R += 100%

function [E, R] = greenhouse_effect_final(modelNum, E_old, dt, iter, prev_Cw, Cw, Cb, Aw, Ab, R)

    % Default constants, if not specified
    alpha = 2;
    beta = 1;

    if modelNum == 0
        delE = 0;

    elseif modelNum == 1
        delE = 1;

    elseif modelNum == 2
        alpha = 5;
        beta = 0.5;
        delE = alpha - beta * E_old * (Aw + Ab);

    elseif modelNum == 3
        delE = alpha * Cw - beta * Cb;

    elseif modelNum == 4
        if iter > 5
            if (Cw - prev_Cw) < 0
                R = R + 0.1;
                if R > 1
                    R = 1;
                end
            else
                R = 0;
            end
        end

        delE = alpha * (1 - R) * Cw - beta * Cb;

    elseif modelNum == 5
        if iter > 5
            if (Cw - prev_Cw) < 0
                R = 1;
            else
                R = 0;
            end
        end

        delE = alpha * (1 - R) * Cw - beta * Cb;

    elseif modelNum == 6
        if iter > 5
            if Cw > (2 * Cb)
                R = R + 0.02;
            else
                R = R - 0.02;
            end
        end

        if R > 1
            R = 1;
        end
        if R < 0
            R = 0;
        end

        delE = alpha * (1 - R) * Cw - beta * Cb;

    elseif modelNum == 7
        if iter > 5
            if Cw > 2 * Cb
                R = R + 0.1;
            else
                R = R - 0.1;
            end
        end

        if R > 1
            R = 1;
        end
        if R < 0
            R = 0;
        end

        delE = alpha * (1 - R) * Cw - beta * Cb;

    elseif modelNum == 8
        if iter > 5
            if Cw > 2 * Cb
                R = 1;
            else
                R = 0;
            end
        end

        delE = alpha * (1 - R) * Cw - beta * Cb;
    end

    E = E_old + dt * delE;

    if modelNum < 4
        R = 0;
    end
end

