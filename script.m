% Load the dataset
data = load('DataSet3.mat');

%Selecting Signal to work with
receivedSignal = data.HighNoise_RxSignal;


%Removing Cyclic Prefixing
cyclicPrefixLength = data.Channel.Length-1;
OFDMSymbolLength = data.OFDM.FFT_Length;
% Total length of each received OFDM symbol (with cyclic prefix)
totalSymbolLength = OFDMSymbolLength + cyclicPrefixLength;

% Reshape the received signal into symbols, assuming each block is an OFDM symbol
numOFDMSymbols = length(receivedSignal) / totalSymbolLength;  % Number of OFDM symbols
% Preallocate for the received symbols without cyclic prefix
rxSymbolsWithoutCP = zeros(int32(data.OFDM.FFT_Length), int32(numOFDMSymbols));

% Remove cyclic prefix for each OFDM symbol
for i = 1:numOFDMSymbols
    % Calculate start and end index for each OFDM symbol
    startIdx = (i - 1) * totalSymbolLength + 1; % Start index for current symbol
    endIdx = startIdx + totalSymbolLength - 1;  % End index for current symbol
    
    % Check to ensure you are within bounds
    if endIdx <= length(receivedSignal)
        ofdmSymbolWithCP = receivedSignal(startIdx:endIdx); % Extract the symbol including CP
        % Remove the cyclic prefix
        rxSymbolsWithoutCP(:, i) = ofdmSymbolWithCP(cyclicPrefixLength + 1:end); % Keep only the part after the CP
    else
        error('Index out of bounds. Check the length of receivedSignal and totalSymbolLength.');
    end
end

%fft
rxSymbolsInFrequencyDomain = fft(rxSymbolsWithoutCP);

%% Channel Estimation With Kalman Filter
%Creating all of the vectors and matrices I will need
% Get the number of columns in rxSymbolsInFrequencyDomain
[numRows, numCols] = size(rxSymbolsInFrequencyDomain);

% Initialize vectors/matrices of size equal to numCols
h = ones(1, numCols);  % Initialize the previous channel estimate (vector)
P = eye(numCols);       % Initialize the covariance matrix as an identity matrix (size numCols x numCols)
xc_matrix = zeros(1,numCols);  % Pre-allocate matrix for efficiency
for i = 1:numCols
    [xc,lags] = xcorr(rxSymbolsWithoutCP(:,i),numRows,'coeff');
    xc_matrix(i) = xc(numRows);      % Store result in the i-th column
end   


rho = mean(abs(xc_matrix(1,:))) .* ones(1, numCols);  % Initialize the correlation coefficient vector (same value for all processes) 0.12 for Dataset1
Q = 0.03 * eye(numCols);     % Initialize process noise covariance (size numCols x numCols)
R = 0.05 * eye(numCols);     % Initialize measurement noise covariance (size numCols x numCols)

% Calculate the time shift effect
N = length(rxSymbolsWithoutCP);
omega = 2 * pi * 1 / N;
phaseShift = exp(-1j * omega);
%phaseShift = exp(-1j * pi/10);

% New measurement data 
x = ones(1, numCols)*data.OFDM.PilotSymbol;  % Transmitted pilot signals
%Storing the filtered symbols
KalmanfilteredSymbols = zeros(numRows,numCols);

for k = 1:numRows
    % Get the current received symbol
    receivedSymbol = rxSymbolsInFrequencyDomain(k, :);  % Current row (subcarrier)

    % Check if the current index matches a pilot signal index
    if ismember(k, data.OFDM.PilotIndices)
        % Update step: Call the update function
        [h, P] = kalman_update(h, P, receivedSymbol, x, R);
    end
    % Prediction step: Call the prediction function
    [h, P] = kalman_prediction(h, P, rho*phaseShift, Q);

    KalmanfilteredSymbols(k,:) = (1./h).*receivedSymbol;
end


%% Channel Estimation Weiener Filter
% Lets use a Wiener Filter of order 3 to estimate the signal
%Creating all of the vectors and matrices I will need
WienerfilteredSymbols = zeros(numRows,numCols);

pilotOnlyVector = zeros(size(rxSymbolsInFrequencyDomain));
pilotOnlyVector(data.OFDM.PilotIndices,:)  = rxSymbolsInFrequencyDomain(data.OFDM.PilotIndices, :);

desiredPilotVector = zeros(size(rxSymbolsInFrequencyDomain));
desiredPilotVector(data.OFDM.PilotIndices,:) = data.OFDM.PilotSymbol;

pilotOnlyVectorInTime = ifft(pilotOnlyVector);
desiredPilotVectorInTime = ifft(desiredPilotVector);

% Weiner filter for time series

for k = 1:numCols 
    % Compute autocorrelation of the received signal
    R_xx = xcorr(pilotOnlyVectorInTime(:,k), 'unbiased'); % Autocorrelation of received signal
    R_xx = R_xx(N:end); % Keep only the second half (positive lags)
    % Compute cross-correlation between desired and received signal
    R_dx = xcorr(desiredPilotVectorInTime(:,k), pilotOnlyVectorInTime(:,k), 'unbiased'); % Cross-correlation
    R_dx = R_dx(N:end); % Keep only the second half (positive lags)

    % Create the autocorrelation matrix
    L = 4; % Length of the Wiener filter
    R_xx_matrix = zeros(L, L);
    
    for i = 1:L
        R_xx_matrix(i, :) = R_xx(i:i + L - 1)'; % Fill the matrix with autocorrelations
    end
    
    %Create the cross-correlation vector
    R_dx_vector = R_dx(1:L); % Take the first L elements for cross-correlation

    h = R_xx_matrix \ R_dx_vector; % Solve for filter coefficients
    H = fft(h,numRows);
    WienerfilteredSymbols(:,k) = H.*rxSymbolsInFrequencyDomain(:,k);
end

%% Choosing What Method to Use
% Uncomment the section you want to use and comment the other 2
%KALMAN FILTER
pilotSymbols = KalmanfilteredSymbols(data.OFDM.PilotIndices, :);
dataSymbols = KalmanfilteredSymbols(data.OFDM.DataIndices, :);

%WEINER FILTER
% pilotSymbols = WienerfilteredSymbols(data.OFDM.PilotIndices, :);
% dataSymbols = WienerfilteredSymbols(data.OFDM.DataIndices, :);

%NO FILTER
% pilotSymbols = rxSymbolsInFrequencyDomain(data.OFDM.PilotIndices, :);
% dataSymbols = rxSymbolsInFrequencyDomain(data.OFDM.DataIndices, :);

%% Plotting

% Plotting pilot symbols received vs transmitted
pilotSymbolsTransmitted = data.OFDM.PilotSymbol * ones(size(pilotSymbols(:,1))); 
figure;
scatter(real(pilotSymbolsTransmitted(:)), imag(pilotSymbolsTransmitted(:)), 'filled');
hold on;
scatter(real(pilotSymbols(:)), imag(pilotSymbols(:)), 'filled');
xlabel('Real Part');
ylabel('Imaginary Part');
title('Transmitted vs Rreceived Pilot Symbols');
legend('Transmitted pilot Symbols', 'Received pilot Symbols');
grid on;
axis equal;
hold off;

% Plotting the difference for every column in PilotSymbols
figure;
hold on;
% Iterate over each column in pilotSymbols
for i = 1:size(pilotSymbols, 2)
    % Calculate the absolute difference for the i-th column
    absDifference = abs(pilotSymbols(:, i) - pilotSymbolsTransmitted);
    
    % Plot the absolute difference for the current column
    plot(data.OFDM.PilotIndices, absDifference(:), 'LineWidth', 1.5, 'DisplayName', ['Column ' num2str(i)]);
end

% Customize the plot
xlabel('Index');
ylabel('Absolute Difference');
title('Absolute Difference Between Received and Transmitted Pilot Symbols');
legend show; % Show legend to identify each column
grid on;
hold off; % Release the hold

%% Decoding
% Define un umbral de tolerancia
tolerance = 1e-6;
nonZeroMask = abs(dataSymbols) > tolerance;
dataSymbolsNoPadding = dataSymbols(nonZeroMask);

rxData = qamdemod(dataSymbolsNoPadding, 4, 'OutputType', 'bit', 'UnitAveragePower', true);
length(dataSymbolsNoPadding)

image = 255*reshape(rxData(1:(data.ImageSize(1)*data.ImageSize(2))),data.ImageSize);
%Plotting the Image
figure;
DecodeImage = cast(image,'uint8');
        imshow(DecodeImage);

