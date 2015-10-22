%% particle filter for 1D signal
clear
N = 500;
iter = 20;
x_s0 = 1;
sigma = 0.5;

X = zeros(1,N);
U = zeros(1,N);
w = zeros(1,N);
X_sample = round(normrnd(x_s0,sigma,1,N)); % 粒子的初始状态
w0 = 1/N * ones(1,N); % 粒子的初始权重

sigma_noise = 0.01;
for i = 1:iter
    noise_state = normrnd(0,sigma_noise,1,N); % 模型噪声，高斯噪声
    X = X_sample + 10 + noise_state; % 状态转移方程 
    noise_measure = normrnd(0,sigma_noise,1,N); % 观测噪声，高斯噪声
    U = X + noise_measure; % 观测方程
    delta = (U - X).^2;
    w = delta/sum(delta);
%     x_new = round(sum(w.*X));
    x_new = sum(w.*X);
    scatter(x_new,i,40,'filled');hold on;axis([0 200 0 25]);
%     X_sample = round(normrnd(x_new,sigma,1,N));
    X_sample = normrnd(x_new,sigma,1,N);
    scatter(X_sample,i*ones(1,N),3,'filled');xlabel('状态值/观测值');ylabel('次数');
end

%% particle filter for 2D signal

    
    





















