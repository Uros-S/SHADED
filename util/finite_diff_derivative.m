function derivative = finite_diff_derivative(signal, t, order)
% This function estimates the derivative of a signal using finite 
% differences. In particular, forward difference is used at the first
% point(s), central difference at interior points, and backward difference
% at the last point(s).

    if order == 4 && length(signal) < 6
        error('Order-4 derivative requires at least 6 points!');
    end

    % Check that inputs have the same size
    if length(signal) ~= length(t)
        error('Signal and t vectors must have the same length!');
    end

    n = length(signal);
    derivative = zeros(size(signal, 1), size(signal, 2));

    if n < 5
        error('Signal must have at least 5 points!');
    end

    % Forward difference at boundary points
    if order == 1
        derivative(1) = (signal(2) - signal(1)) / (t(2) - t(1));

    elseif order == 2
        derivative(1) = (signal(3) - 2*signal(2) + signal(1)) / (t(2) - t(1))^2;

    elseif order == 3
        derivative(1) = (signal(4) - 3*signal(3) + 3*signal(2) - signal(1)) / (t(2) - t(1))^3;
        derivative(2) = (signal(5) - 3*signal(4) + 3*signal(3) - signal(2)) / (t(3) - t(2))^3;

    elseif order == 4
        % Two forward-difference boundary points (need indices 1..5 and 2..6)
        derivative(1) = (signal(1) - 4*signal(2) + 6*signal(3) - 4*signal(4) + signal(5)) ...
                        / (t(2) - t(1))^4;
        derivative(2) = (signal(2) - 4*signal(3) + 6*signal(4) - 4*signal(5) + signal(6)) ...
                        / (t(3) - t(2))^4;
    end

    % Central difference at interior points 
    for i = 2:n-1
        if order == 1
            derivative(i) = (signal(i+1) - signal(i-1)) / (t(i+1) - t(i-1));

        elseif order == 2
            derivative(i) = (signal(i+1) - 2*signal(i) + signal(i-1)) / (t(i) - t(i-1))^2;

        elseif order == 3
            if i == 2 || i == n-1
                continue;  % filled by forward/backward boundary blocks
            end
            derivative(i) = (signal(i+2) - 2*signal(i+1) + 2*signal(i-1) - signal(i-2)) ...
                            / (2*(t(i) - t(i-1))^3);

        elseif order == 4
            if i <= 2 || i >= n-1
                continue;  % filled by forward/backward boundary blocks
            end
            derivative(i) = (signal(i+2) - 4*signal(i+1) + 6*signal(i) - 4*signal(i-1) + signal(i-2)) ...
                            / (t(i) - t(i-1))^4;
        end
    end

    % Backward difference at boundary points
    if order == 1
        derivative(n) = (signal(n) - signal(n-1)) / (t(n) - t(n-1));

    elseif order == 2
        derivative(n) = (signal(n) - 2*signal(n-1) + signal(n-2)) / (t(n) - t(n-1))^2;

    elseif order == 3
        derivative(n-1) = (signal(n-1) - 3*signal(n-2) + 3*signal(n-3) - signal(n-4)) ...
                          / (t(n-1) - t(n-2))^3;
        derivative(n)   = (signal(n)   - 3*signal(n-1) + 3*signal(n-2) - signal(n-3)) ...
                          / (t(n) - t(n-1))^3;

    elseif order == 4
        % Two backward-difference boundary points
        derivative(n-1) = (signal(n-1) - 4*signal(n-2) + 6*signal(n-3) - 4*signal(n-4) + signal(n-5)) ...
                          / (t(n-1) - t(n-2))^4;
        derivative(n)   = (signal(n)   - 4*signal(n-1) + 6*signal(n-2) - 4*signal(n-3) + signal(n-4)) ...
                          / (t(n) - t(n-1))^4;
    end

end