% kalman_prediction: Prediction step of the Kalman filter
    %
    % Inputs:
    % h_prev : Nx1 vector of previous channel estimates
    % P_prev : NxN previous error covariance matrix
    % p      : Nx1 vector of correlation coefficients
    % Q      : NxN process noise covariance matrix
    %
    % Outputs:
    % h_pred : Nx1 predicted channel estimate
    % P_pred : NxN predicted error covariance matrix
function [h_pred, P_pred] = kalman_prediction(h_prev, P_prev, p, Q)
    % Predicted state
    h_pred = p .* h_prev;  % Element-wise multiplication for each process

    % Predicted covariance
    P_pred = diag(p.^2) * P_prev + Q;  % Include process noise
end
