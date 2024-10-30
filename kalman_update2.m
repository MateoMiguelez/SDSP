% kalman_update: Update step of the Kalman filter
    %
    % Inputs:
    % h_pred : Nx1 predicted channel estimate
    % P_pred : NxN predicted error covariance matrix
    % z      : Nx1 vector of received signals (observations)
    % x      : Nx1 vector of pilot signals
    % R      : NxN measurement noise covariance matrix
    %
    % Outputs:
    % h_upd  : Nx1 updated channel estimate
    % P_upd  : NxN updated error covariance matrix

function [h_upd, P_upd] = kalman_update2(h_pred, P_pred, z, x, R)
    
    % Measurement innovation (difference between predicted and received signal)
    y_k = z - (h_pred .* x);  % Innovation

    % Innovation covariance
    S = diag(P_pred) .* (x.^2) + R;  % Element-wise operation with x and diagonal of P_pred

    % Kalman gain
    K = (diag(P_pred) .* x) ./ S;  % Element-wise division

    % Updated state estimate
    h_upd = h_pred + y_k*K;  % Update estimate

    % Updated covariance estimate
    P_upd = (eye(size(P_pred)) - diag(K .* x)) * P_pred;  % Update covariance
end
