clc; clear
%% ----------------------------
% 1. Parameters
%% ----------------------------
N = 8;          % Polar block length
K = 4;          % Number of information bits
SNR_dB = 5;     % Channel SNR in dB (the bigger the better)
A = 2;          % Amplitude (increase to reduce error)

frozen_idx = [1:3 5];  % First 4 bits frozen
info_idx = [4 6:8];    % Last 4 bits are info bits

message = 'hello';

%% ----------------------------
% 2. Convert message to binary
%% ----------------------------
bin_message = '';
for i = 1:length(message)
    bin_message = [bin_message, dec2bin(double(message(i)),8)];
end
disp(['Binary message: ', bin_message]);

%% ----------------------------
% 3. Generate Polar Transform Matrix (Kronecker)
%% ----------------------------
F = [1 0; 1 1];
F_N = 1;
n = log2(N);
for i = 1:n
    F_N = kron(F_N,F);
end
F_N = mod(F_N,2);
%% ----------------------------
% 4. Divide message into K-bit blocks
%% ----------------------------

num_blocks = length(bin_message)/K;

disp(['Number of blocks: ', num2str(num_blocks)]);
%% ----------------------------
% 5. Simulation loop
%% ----------------------------
decoded_bits = [];
sigma = sqrt(1/(2*10^(SNR_dB/10))); % Noise std for AWGN
for blk = 1:num_blocks
    % Extract info bits for this block
    info_bits = bin_message((blk-1)*K+1 : blk*K);
    % Build N-bit vector with frozen bits = 0
    u = zeros(1,N);
    u(info_idx) = info_bits - '0';  % Convert char to 0/1
    % Polar encoding
    x = mod(u * F_N,2);
    % BPSK modulation
    tx = A*(1 - 2*x); % 0->+1, 1->-1
    % AWGN channel
    y = tx + sigma*randn(1,N);
    % LLR computation
    LLR = 2*y/(sigma^2);
    % Simple SC decoding (hard LLR)
    u_hat = zeros(1,N);
    for i = 1:N
        if ismember(i,frozen_idx)
            u_hat(i) = 0;
        else
            u_hat(i) = LLR(i)<0; % 0 if LLR>0, 1 if LLR<0
        end
    end
    for i= 1:N-1
      if ismember(i,info_idx)
        u_hat(i) = xor(u_hat(i),u_hat(8));
      end
    end
    % Collect decoded info bits
    decoded_bits = [decoded_bits u_hat(info_idx)];
end
%% ----------------------------
% 6. Convert decoded bits to string
%% ----------------------------
decoded_message = '';
for i = 1:8:length(decoded_bits)
    char_bits = decoded_bits(i:i+7);
    decoded_message = [decoded_message char(bin2dec(num2str(char_bits)))];
end
disp(['Decoded message: ', decoded_message]);

