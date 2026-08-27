fprintf('\n===== A_METHOD_C0_VGJ_RANGE_V11_FINAL_RUNNING =====\n');
clc;
clearvars;
fprintf('\n===== V8_RELERR_1003_OVERALL_WINDOW_ONLY_RUNNING =====\n');
clear functions;
close all force;
fclose('all');

% clc;
% clear;
% close all;
% 
% % 加载数据
% data = importdata('Final_data.txt');
% 
% % 读取完整数据
% data1 = importdata('puredata2.txt');
% data1 = data1.data;
% 
% % 找出 data1 中哪些行不在 data 中
% % 首先将每一行转为字符串形式以便比较
% str_data1 = string(num2str(data1));
% str_data = string(num2str(data));
% 
% % 找出 data1 中每一行是否出现在 data 中
% [~, ia] = ismember(str_data1, str_data);
% 
% % ia==0 表示不在 data 中的索引（即不同的数据）
% diff_idx = find(ia == 0);
% 
% % 显示不同数据的索引
% disp('data1 中与 data 不同的行索引是：');
% disp(diff_idx);
% 
% % 可选：查看这些不同的数据
% different_data = data1(diff_idx, :);
clc; clear; clear functions; close all;
fprintf('\n===== A_METHOD_C0_VGJ_RANGE_V10_CLEAR_OUTPUT_RUNNING =====\n');
% 导入数据
data = importdata('Final_data.txt');
% 基本参数
rho_f = data(:,11);    % 液体密度 (kg/m^3)
rho_g = data(:,12);    % 气体密度 (kg/m^3)
j_f = data(:,1);       % 液相速度 (m/s)
j_g = data(:,2);       % 气相速度 (m/s)
d = data(:,8);       % 直径 (m)
p = data(:,9);       % pitch (m)
mu_f = data(:,13);       % 液体动力粘度 (Pa·s)
mu_g = data(:,14);       % 气体动力粘度 (Pa·s)
Re_f = data(:,15);
Re_g = data(:,16);
alpha_m = data(:,6);   % 实测空泡份额 (-)
j_total = j_f + j_g;
g = 9.81;              % 重力加速度 (m/s^2)
P=data(:,5)*1e6;          %系统压力 (Pa)
G=data(:,17);
N = size(data,1);
% 使用 REFPROP 获取表面张力（单位 N/m）
T = data(:,4)+ 273.15;  % 温度 (K)
sigma = zeros(N,1);
for i = 1:N
    [~, sigma(i)] = refpropm('SIGMA', 'T', T(i), 'Q', 0, 'Water');  % 提取正确的 σ
end
Dh =data(:,10);         % 水力直径 (m)
% 公共变量
sqrt_ratio = sqrt(rho_f ./ rho_g);

%% Jowitt et al 模型
C0_Jowitt = 1 + 0.796 * exp(-0.061 * sqrt_ratio);
v_gj_Jowitt = 0.034 * (sqrt_ratio - 1);
alpha_Jowitt = j_g ./ (C0_Jowitt .* j_total + v_gj_Jowitt);
m_rel_Jowitt = mean(abs((alpha_Jowitt - alpha_m) ./ alpha_m));

%% Morooka 模型
N = length(j_total);  % 总数据点数
C0_Morooka = ones(N,1);
v_gj_Morooka = zeros(N,1);

% Bond 数计算
Bo = g .* sigma .* (rho_f - rho_g) ./ rho_f.^2;

% 高速 / 低速工况划分
high_flow = j_total >= 5;
low_flow  = ~high_flow;

% 高速工况
C0_Morooka(high_flow) = 1.08;
v_gj_Morooka(high_flow) = 3.04 * (Bo(high_flow)).^0.25;

% 低速工况
C0_Morooka(low_flow) = 1.13;
v_gj_Morooka(low_flow) = 1.41 * (Bo(low_flow)).^0.25;

% 计算 void fraction
alpha_Morooka = j_g ./ (C0_Morooka .* j_total + v_gj_Morooka);

% 误差评估
rel_err_high = mean(abs((alpha_Morooka(high_flow) - alpha_m(high_flow)) ./ alpha_m(high_flow)));
rel_err_low  = mean(abs((alpha_Morooka(low_flow) - alpha_m(low_flow)) ./ alpha_m(low_flow)));

% 加权平均相对误差
N_high = sum(high_flow);
N_low = sum(low_flow);
m_rel_Morooka = (rel_err_high * N_high + rel_err_low * N_low) / (N_high + N_low);

%% Bestion 模型 - 保存所有 vgj
counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
start_idx = 1;
total_N = sum(counts);
% 储存每段误差
segment_errors = zeros(length(counts), 1);
% 初始化整体 vgj 向量
v_gj_Bestion_all = zeros(total_N, 1);
alpha_Bestion_all = zeros(total_N, 1);

for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;
    % 当前段数据
    rho_f_seg = rho_f(idx_range);
    rho_g_seg = rho_g(idx_range);
    j_g_seg   = j_g(idx_range);
    j_total_seg = j_total(idx_range);
    alpha_m_seg = alpha_m(idx_range);
    Dh_seg    = Dh(idx_range);
    % Bestion 模型计算
    C0_Bestion      = 1.2 - 0.2 * sqrt(rho_g_seg ./ rho_f_seg);
    v_gj_Bestion    = 0.188 * sqrt(g .* Dh_seg .* (rho_f_seg - rho_g_seg) ./ rho_g_seg);
    alpha_Bestion   = j_g_seg ./ (C0_Bestion .* j_total_seg + v_gj_Bestion);
    % 误差
    segment_errors(k) = mean(abs((alpha_Bestion - alpha_m_seg) ./ alpha_m_seg), 'omitnan');
    % 保存所有结果
    v_gj_Bestion_all(idx_range) = v_gj_Bestion;
    alpha_Bestion_all(idx_range) = alpha_Bestion;
    start_idx = start_idx + counts(k);
end
% 加权平均误差
m_rel_Bestion = sum(segment_errors .* counts') / total_N;


%% Chexal 模型 
counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
start_idx = 1;
total_N = sum(counts);

% 初始化全体保存变量
alpha_chexal_all = NaN(total_N,1);
vgj_chexal_all   = NaN(total_N,1);
C0_chexal_all    = NaN(total_N,1);
segment_errors_chexal = NaN(length(counts),1);

max_iter = 400;
tol = 1e-8;
P_crit = 22.064e6;  % 水的临界压力 (Pa)
for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;
    for ii = 1:length(idx_range)
        i = idx_range(ii);
        % === 基本输入 ===
        rho_fi = rho_f(i);
        rho_gi = rho_g(i);
        mu_fi  = mu_f(i);
        mu_gi  = mu_g(i);  % 可省略
        Dhi    = Dh(i);
        sigma_i = sigma(i);
        Pi = P(i);
        jgi = j_g(i);
        jfi = j_f(i);
        jtoti = jfi + jgi;
        Re_fi = Re_f(i);
        Re_gi = Re_g(i);
        Rei = max(Re_fi, Re_gi);
        alpha_old = max(alpha_m(i), 0.01);
        for iter = 1:max_iter
            % ==== 模型参数 ====
            C1 = 4 * P_crit^2 / (P_crit - Pi);
            C5 = sqrt(150 * rho_gi / rho_fi);
            C7 = (0.09144 / Dhi)^0.6;
            % C2
            rho_ratio = rho_fi / rho_gi;
            if rho_ratio <= 18
                C2 = 0.4757 * (log(rho_ratio))^0.7;
            else
                if C5 >= 1
                    C2 = 1;
                else
                    C2 = (1 - exp(-C5 / (1 - C5))).^-1;
                end
            end
            
            % C3
            C3 = max(0.5, 2 * exp(-Re_fi / 60000));
            
            % C4
            if C7 >= 1
                C4 = 1;
            else
                C4 = 1 / (1 - exp(-C7 / (1 - C7)));
            end
            % B1, r, K0
            B1 = min(0.8, 1 / (1 + exp(-Rei / 60000)));
            r = (1 + 1.57 * (rho_gi / rho_fi)) / (1 - B1);
            K0 = B1 + (1 - B1) * (rho_gi / rho_fi)^0.25;
            % C9
            C9 = (1 - alpha_old)^B1;

            % L(alpha) 蒸汽-水混合 - 分段使用不同公式
            if i <= 343
                L = (1 - exp(-C1 * alpha_old)) / (1 - exp(-C1));
            else
                L = min(1.15 * alpha_old^0.45, 1);
            end

            % C0 & vgj
            C0 = L / (K0 + (1 - K0) * alpha_old^r);
            vgj = 1.41 * (sigma_i * g * (rho_fi - rho_gi) / rho_fi^2)^0.25 * C2 * C3 * C4 * C9;
            % alpha 更新
            denom = C0 * jtoti + vgj;
            if denom <= 0
                break;
            end
            alpha_new = jgi / denom;
            if ~isfinite(alpha_new) || alpha_new < 0 || alpha_new > 1
                break;
            end
            
            if abs(alpha_new - alpha_old) < tol
                break;
            end
            alpha_old = alpha_new;
        end
        
        % 保存结果
        alpha_chexal_all(i) = alpha_new;
        vgj_chexal_all(i) = vgj;
        C0_chexal_all(i) = C0;
    end
    
    % 当前段误差
    alpha_m_seg = alpha_m(idx_range);
    valid_idx = alpha_m_seg > 0.001 & isfinite(alpha_chexal_all(idx_range));
    rel_err = abs((alpha_chexal_all(idx_range(valid_idx)) - alpha_m_seg(valid_idx)) ./ alpha_m_seg(valid_idx));
    segment_errors_chexal(k) = mean(rel_err, 'omitnan');
    
    start_idx = start_idx + counts(k);
end


% 总体加权误差
m_rel_chexal = sum(segment_errors_chexal .* counts') / total_N;
v_gj_chexal=vgj_chexal_all;

%% Maier and Coddington 模型
C0_MC = 2.57e-3 .* P*1e-6 + 1.0062;

v_gj_MC = (6.73e-7 .* (P*1e-6).^2 - 8.81e-5 .* (P*1e-6) + 1.05e-3) .* G ...
    + (5.63e-3 .* (P*1e-6).^2 - 1.23e-1 .* (P*1e-6) + 0.8);

alpha_MC = j_g ./ (C0_MC .* j_total + v_gj_MC);

% 分别计算高速/低速区的相对误差
m_rel_MC = mean(abs((alpha_MC - alpha_m) ./ alpha_m));


%% Julia et al. 模型迭代版（分段算法）

counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
start_idx = 1;

segment_errors_Julia     = zeros(length(counts), 1);
final_iter_error_Julia   = zeros(length(counts), 1);
iteration_counts_Julia   = zeros(length(counts), 1);

alpha_all_model = [];
alpha_all_measured = [];

max_iter = 100;
tol = 1e-6;

for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;
    N_seg = length(idx_range);
    
    % 分段数据
    jg_seg = j_g(idx_range);
    jt_seg = j_total(idx_range);
    alpha_m_seg = alpha_m(idx_range);
    rho_f_seg = rho_f(idx_range);
    rho_g_seg = rho_g(idx_range);
    sigma_seg = sigma(idx_range);
    d_seg = d(idx_range);
    p_seg = p(idx_range);
    
    d_over_p = d_seg ./ p_seg;
    rho_ratio = rho_g_seg ./ rho_f_seg;
    sqrt_term = sqrt(2) * (sigma_seg .* g .* (rho_f_seg - rho_g_seg) ./ rho_f_seg.^2).^(0.25);
    
    alpha_seg = NaN(N_seg,1);
    final_diff_seg = NaN(N_seg,1);
    
    for i = 1:N_seg
        jg_i = jg_seg(i);
        jt_i = jt_seg(i);
        alpha_init = alpha_m_seg(i);
        if alpha_init < 0.01 || alpha_init > 0.9 || isnan(alpha_init)
            alpha_old = 0.2;
        else
            alpha_old = alpha_init;
        end
        
        % 基本量
        ro_ratio_sqrt = sqrt(rho_ratio(i));
        d_by_p = min(max(d_over_p(i), 0.3), 0.7);
        d_i = d_seg(i);
        p_i = p_seg(i);
        sigma_i = sigma_seg(i);
        rf_i = rho_f_seg(i);
        rg_i = rho_g_seg(i);
        
        % === Lmax, Db, Bsf ===
        Lmax = sqrt(2) * p_i - d_i;
        if Lmax <= 0 || isnan(Lmax)
            continue;
        end
        Db = sqrt(0.4 * sigma_i / (g * max(rf_i - rg_i, 1e-6)));
        D_ratio = Db / Lmax;
        if D_ratio < 0.6
            Bsf = 1 - Db / (0.9 * Lmax);
        else
            Bsf = 0.120 / (D_ratio^2);
        end
        sqrt_prefactor = sqrt_term(i);
        
        % === 迭代 ===
        for iter = 1:max_iter
            C0_vals = [
                (1.03 - 0.03 * ro_ratio_sqrt) * (1 - exp(-26.3 * alpha_old^0.780));
                (1.04 - 0.04 * ro_ratio_sqrt) * (1 - exp(-21.2 * alpha_old^0.762));
                (1.05 - 0.05 * ro_ratio_sqrt) * (1 - exp(-34.1 * alpha_old^0.925));
                ];
            C0 = interp1([0.3;0.5;0.7], C0_vals, d_by_p);
            
            if ~isfinite(C0) || C0 <= 0 || alpha_old >= 1
                break;
            end
            
            vgj = Bsf * sqrt_prefactor * (1 - alpha_old)^1.75;
            denom = C0 * jt_i + vgj;
            if denom <= 0
                break;
            end
            alpha_new = jg_i / denom;
            
            if ~isfinite(alpha_new) || alpha_new < 0 || alpha_new > 1
                break;
            end
            
            diff = abs(alpha_new - alpha_old);
            if diff < tol
                final_diff_seg(i) = diff;
                iteration_counts_Julia(k) = max(iteration_counts_Julia(k), iter);
                break;
            end
            if iter == max_iter
                final_diff_seg(i) = diff;
                iteration_counts_Julia(k) = max_iter;
            end
            alpha_old = alpha_new;
        end
        alpha_seg(i) = alpha_new;
    end
    
    % === 每段计算 ===
    valid_idx = alpha_m_seg > 1e-3 & isfinite(alpha_seg) & alpha_seg >= 0 & alpha_seg <= 1;
    rel_err_seg = abs((alpha_seg(valid_idx) - alpha_m_seg(valid_idx)) ./ alpha_m_seg(valid_idx));
    
    segment_errors_Julia(k) = mean(rel_err_seg, 'omitnan');
    final_iter_error_Julia(k) = mean(final_diff_seg(valid_idx), 'omitnan');
    
    alpha_all_model = [alpha_all_model; alpha_seg];
    alpha_all_measured = [alpha_all_measured; alpha_m_seg];
    
    start_idx = start_idx + counts(k);
end

% === 总体计算 ===
valid_all = alpha_all_measured > 1e-3 & isfinite(alpha_all_model) & alpha_all_model >= 0 & alpha_all_model <= 1;
m_rel_Julia = mean(abs((alpha_all_model(valid_all) - alpha_all_measured(valid_all)) ./ alpha_all_measured(valid_all)));

% ========================================================================
%  关键修正：立即保存 Julia et al. 的 alpha 预测值
%  原因：后面的 Clark 模型也会使用变量名 alpha_all_model，
%        如果这里不保存，Julia 结果会被 Clark 覆盖，导致 ML 收集阶段
%        出现“Julia et al. 未找到对应 alpha 变量”。
% ========================================================================
alpha_Julia_saved = alpha_all_model(:);
alpha_julia_all   = alpha_Julia_saved;
alpha_Julia       = alpha_Julia_saved;

fprintf('Julia et al. alpha 已保存为 alpha_Julia_saved，可参与 PIE-DF / PIE-DF-E。\n');

%% Paranjape et al. 模型
C0_Paranjape = 1.05;
v_gj_Paranjape = 0.123;  % m/s, 常数

alpha_Paranjape = j_g ./ (C0_Paranjape .* j_total + v_gj_Paranjape);
m_rel_Paranjape = mean(abs((alpha_Paranjape - alpha_m) ./ alpha_m));


%% Kamei et al. 模型
C0_Kamei = ones(N,1);
v_gj_Kamei = ones(N,1);

low_jg = j_g < 1.5;
high_jg = ~low_jg;

% 分段设置
C0_Kamei(low_jg) = 1.44;
v_gj_Kamei(low_jg) = 0.58;

C0_Kamei(high_jg) = 0.96;
v_gj_Kamei(high_jg) = 1.32;

% 计算预测值和误差
alpha_Kamei = j_g ./ (C0_Kamei .* j_total + v_gj_Kamei);
m_rel_Kamei = mean(abs((alpha_Kamei - alpha_m) ./ alpha_m));

%% Chen 模型（迭代 + 分段，复用 Clark 的 vgj 结构）
counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
start_idx = 1;
total_N = sum(counts);

alpha_chen_all = NaN(total_N,1);
C0_chen_all    = NaN(total_N,1);
vgj_chen_all   = NaN(total_N,1);
segment_errors_chen = NaN(length(counts),1);

max_iter = 400;
tol = 1e-8;

% 非维量
cap = (sigma .* g .* (rho_f - rho_g) ./ rho_f.^2).^(0.25);
j_g_plus = j_g ./ cap;
j_f_plus = j_f ./ cap;
sqrt_term = sqrt(2) * cap;

for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;
    for ii = 1:length(idx_range)
        i = idx_range(ii);
        
        % 参数提取
        rho_fi = rho_f(i);
        rho_gi = rho_g(i);
        mu_fi = mu_f(i);
        jgi = j_g(i);
        jfi = j_f(i);
        jtoti = jgi + jfi;
        Dhi = Dh(i);
        sigma_i = sigma(i);
        jg_plus_i = j_g_plus(i);
        cap_i = cap(i);
        delta_rho = rho_fi - rho_gi;
        
        L_star = Dhi / sqrt(sigma_i / (g * delta_rho));
        N_mu = mu_fi / sqrt(sigma_i * rho_fi * sqrt(sigma_i / (g * delta_rho)));
        r_term = (rho_gi / rho_fi)^(-0.157);
        
        alpha_old = max(alpha_m(i), 0.0001);
        if alpha_old > 0.99 || isnan(alpha_old)
            alpha_old = 0.0001;
        end
        
        for iter = 1:max_iter
            % === Chen 模型的 C_inf 与 C0 ===
            if jg_plus_i <= 0.5
                C_inf = 4.79 * jg_plus_i + 1;
            else
                C_inf = 3.45 * exp(-0.52 * jg_plus_i^0.51) + 1;
            end
            C0 = C_inf - (C_inf - 1) * sqrt(rho_gi / rho_fi);
            
            % === Clark 模型的 vgj 结构（直接复用）===
            if N_mu <= 0.00225 && L_star <= 30
                vgj_P = 0.0019 * L_star^0.809 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu <= 0.00225 && L_star > 30
                vgj_P = 0.030 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu > 0.00225 && L_star > 30
                vgj_P = 0.92 * r_term * cap_i;
            else
                vgj_P = NaN;
            end
            
            vgj_B = sqrt_term(i) * (1 - alpha_old)^1.75;
            vgj = vgj_B * exp(-1.39 * jg_plus_i) + vgj_P * (1 - exp(-1.39 * jg_plus_i));
            
            % === alpha 更新 ===
            denom = C0 * jtoti + vgj;
            if denom <= 0
                break;
            end
            
            alpha_new = jgi / denom;
            if ~isfinite(alpha_new) || alpha_new < 0 || alpha_new > 1
                break;
            end
            
            if abs(alpha_new - alpha_old) < tol
                break;
            end
            
            alpha_old = alpha_new;
        end
        
        % 保存结果
        alpha_chen_all(i) = alpha_new;
        C0_chen_all(i) = C0;
        vgj_chen_all(i) = vgj;
    end
    
    % 每段误差
    alpha_seg = alpha_chen_all(idx_range);
    alpha_m_seg = alpha_m(idx_range);
    valid_idx = alpha_m_seg > 1e-3 & isfinite(alpha_seg);
    rel_err = abs((alpha_seg(valid_idx) - alpha_m_seg(valid_idx)) ./ alpha_m_seg(valid_idx));
    segment_errors_chen(k) = mean(rel_err, 'omitnan');
    
    start_idx = start_idx + counts(k);
end

% 总加权平均误差
m_rel_Chen = sum(segment_errors_chen .* counts') / total_N;
v_gj_Chen=vgj_chen_all;
C0_chen=C0_chen_all;


%% Ozaki & Hibiki 模型（No.1）- 分段迭代实现
counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
start_idx = 1;
total_N = sum(counts);

alpha_ozaki_all = NaN(total_N,1);
vgj_ozaki_all   = NaN(total_N,1);
C0_ozaki_all    = NaN(total_N,1);
segment_errors_ozaki = NaN(length(counts),1);

max_iter = 400;
tol = 1e-8;

% 非维量
cap = (sigma .* g .* (rho_f - rho_g) ./ rho_f.^2).^(0.25);
j_g_plus = j_g ./ cap;
sqrt_term = sqrt(2) * cap;

for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;
    for ii = 1:length(idx_range)
        i = idx_range(ii);
        
        % 提取参数
        rho_fi = rho_f(i);
        rho_gi = rho_g(i);
        mu_fi  = mu_f(i);
        jgi = j_g(i);
        jfi = j_f(i);
        jtoti = jgi + jfi;
        Dhi = Dh(i);
        sigma_i = sigma(i);
        jg_plus_i = j_g_plus(i);
        cap_i = cap(i);
        delta_rho = rho_fi - rho_gi;
        
        L_star = Dhi / sqrt(sigma_i / (g * delta_rho));
        N_mu = mu_fi / sqrt(sigma_i * rho_fi * sqrt(sigma_i / (g * delta_rho)));
        r_term = (rho_gi / rho_fi)^(-0.157);
        
        % 初始 alpha
        alpha_old = max(alpha_m(i), 0.0001);
        if alpha_old > 0.99 || isnan(alpha_old)
            alpha_old = 0.0001;
        end
        
        for iter = 1:max_iter
            % === Ozaki 模型的 C0 ===
            C0 = 1.08 - 0.08 * sqrt(rho_gi / rho_fi);
            C0_ozaki_all(i) = C0;
            
            % === Clark vgj 结构复用 ===
            if N_mu <= 0.00225 && L_star <= 30
                vgj_P = 0.0019 * L_star^0.809 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu <= 0.00225 && L_star > 30
                vgj_P = 0.030 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu > 0.00225 && L_star > 30
                vgj_P = 0.92 * r_term * cap_i;
            else
                vgj_P = NaN;
            end
            
            vgj_B = sqrt_term(i) * (1 - alpha_old)^1.75;
            vgj = vgj_B * exp(-1.39 * jg_plus_i) + vgj_P * (1 - exp(-1.39 * jg_plus_i));
            vgj_ozaki_all(i) = vgj;
            
            denom = C0 * jtoti + vgj;
            if denom <= 0
                break;
            end
            
            alpha_new = jgi / denom;
            if ~isfinite(alpha_new) || alpha_new < 0 || alpha_new > 1
                break;
            end
            
            if abs(alpha_new - alpha_old) < tol
                break;
            end
            
            alpha_old = alpha_new;
        end
        
        alpha_ozaki_all(i) = alpha_new;
    end
    
    % 分段误差
    alpha_seg = alpha_ozaki_all(idx_range);
    alpha_m_seg = alpha_m(idx_range);
    valid_idx = alpha_m_seg > 1e-3 & isfinite(alpha_seg);
    rel_err = abs((alpha_seg(valid_idx) - alpha_m_seg(valid_idx)) ./ alpha_m_seg(valid_idx));
    segment_errors_ozaki(k) = mean(rel_err, 'omitnan');
    
    start_idx = start_idx + counts(k);
end

% 加权平均误差
m_rel_Ozaki = sum(segment_errors_ozaki .* counts') / total_N;
v_gj_0zaki = vgj_ozaki_all;



%% Clark et al. 模型（2014，迭代 + 分段修复版）
counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
start_idx = 1;

segment_errors_Clark     = zeros(length(counts), 1);
final_iter_error_Clark   = zeros(length(counts), 1);
iteration_counts_Clark   = zeros(length(counts), 1);

alpha_all_model = [];
alpha_all_measured = [];

max_iter = 400;
tol = 1e-8;
C0_all = [];

for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;
    N_seg = length(idx_range);
    
    % === 分段数据 ===
    jg_seg = j_g(idx_range);
    jf_seg = j_f(idx_range);
    jt_seg = j_total(idx_range);
    alpha_m_seg = alpha_m(idx_range);
    rho_f_seg = rho_f(idx_range);
    rho_g_seg = rho_g(idx_range);
    mu_f_seg = mu_f(idx_range);
    sigma_seg = sigma(idx_range);
    Dh_seg = Dh(idx_range);
    
    delta_rho_seg = rho_f_seg - rho_g_seg;
    cap_seg = (sigma_seg .* g .* delta_rho_seg ./ rho_f_seg.^2).^(0.25);
    jg_plus = jg_seg ./ cap_seg;
    jf_plus = jf_seg ./ cap_seg;
    j_total_plus = jt_seg ./ cap_seg;
    
    sqrt_term_seg = sqrt(2) * cap_seg;
    
    alpha_seg = NaN(N_seg,1);
    final_diff_seg = NaN(N_seg,1);
    
    C0_seg = NaN(N_seg, 1);
    
    for i = 1:N_seg
        alpha_old = max(alpha_m_seg(i), 0.0001);
        if alpha_old > 0.99 || isnan(alpha_old)
            alpha_old = 0.0001;
        end
        
        C0 = 1;  % 初始猜测
        
        for iter = 1:max_iter
            D_star = Dh_seg(i) / sqrt(sigma_seg(i) / (g * delta_rho_seg(i)));
            N_mu   = mu_f_seg(i) / sqrt(sigma_seg(i) * rho_f_seg(i) * sqrt(sigma_seg(i) / (g * delta_rho_seg(i))));
            r      = (rho_g_seg(i) / rho_f_seg(i))^(-0.157);
            
            if N_mu <= 0.00225 && D_star <= 30
                vgj_P = 0.0019 * D_star^0.809 * r * N_mu^(-0.562)*cap_seg(i);
            elseif N_mu <= 0.00225 && D_star > 30
                vgj_P = 0.030 * r * N_mu^(-0.562)*cap_seg(i);
            elseif N_mu > 0.00225 && D_star > 30
                vgj_P = 0.92 * r*cap_seg(i);
            else
                vgj_P = NaN;
            end
            
            vgj_B    = sqrt_term_seg(i) * (1 - alpha_old)^1.75;
            vgj      = vgj_B * exp(-1.39 * jg_plus(i)) + vgj_P * (1 - exp(-1.39 * jg_plus(i)));
            vgj_plus = vgj / cap_seg(i);
            
            alpha_crit = min(0.0284 * jf_plus(i) + 0.125, 0.52);
            m = 1 / (1 - alpha_crit * C0);
            b = vgj_plus * alpha_crit / (1 - alpha_crit * C0);
            jplus_coo_max = m * jf_plus(i) + b;
            
            C_inf_H_current = 1.1 + 1.84 * exp(-0.100 * j_total_plus(i));
            C_inf_H_coomax  = 1.1 + 1.84 * exp(-0.100 * jplus_coo_max);
            
            if j_total_plus(i) <= jplus_coo_max
                C_inf = ((C_inf_H_coomax - 1) / (jplus_coo_max - jf_plus(i))) * jg_plus(i) + 1;
            else
                C_inf = C_inf_H_current;
            end
            
            C0 = C_inf - (C_inf - 1) * sqrt(rho_g_seg(i) ./ rho_f_seg(i));
            C0_seg(i) = C0;
            
            denom_alpha = C0 * jt_seg(i) + vgj;
            if denom_alpha <= 0
                break;
            end
            
            alpha_new = jg_seg(i) / denom_alpha;
            if ~isfinite(alpha_new) || alpha_new < 0 || alpha_new > 1
                break;
            end
            
            diff = abs(alpha_new - alpha_old);
            if diff < tol
                final_diff_seg(i) = diff;
                iteration_counts_Clark(k) = max(iteration_counts_Clark(k), iter);
                break;
            end
            if iter == max_iter
                final_diff_seg(i) = diff;
                iteration_counts_Clark(k) = max_iter;
            end
            
            alpha_old = alpha_new;
        end
        
        alpha_seg(i) = alpha_new;
    end
    
    valid_idx = alpha_m_seg > 1e-3 & isfinite(alpha_seg) & alpha_seg >= 0 & alpha_seg <= 1;
    rel_err_seg = abs((alpha_seg(valid_idx) - alpha_m_seg(valid_idx)) ./ alpha_m_seg(valid_idx));
    segment_errors_Clark(k) = mean(rel_err_seg, 'omitnan');
    final_iter_error_Clark(k) = mean(final_diff_seg(valid_idx), 'omitnan');
    
    alpha_all_model = [alpha_all_model; alpha_seg];
    alpha_all_measured = [alpha_all_measured; alpha_m_seg];
    C0_all = [C0_all; C0_seg];
    
    start_idx = start_idx + counts(k);
end

m_rel_Clark = sum(segment_errors_Clark .* counts') / sum(counts);

%% Ozaki and Hibiki 模型 (No.1 分段加权)
counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
start_idx = 1;
total_N = sum(counts);

alpha_ozaki_all = NaN(total_N,1);
vgj_ozaki_all   = NaN(total_N,1);
C0_ozaki_all    = NaN(total_N,1);
segment_errors_ozaki = NaN(length(counts),1);

max_iter = 400;
tol = 1e-8;

% 非维量
cap = (sigma .* g .* (rho_f - rho_g) ./ rho_f.^2).^(0.25);
j_g_plus = j_g ./ cap;
sqrt_term = sqrt(2) * cap;

for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;
    for ii = 1:length(idx_range)
        i = idx_range(ii);
        
        % 提取参数
        rho_fi = rho_f(i);
        rho_gi = rho_g(i);
        mu_fi  = mu_f(i);
        jgi = j_g(i);
        jfi = j_f(i);
        jtoti = jgi + jfi;
        Dhi = Dh(i);
        sigma_i = sigma(i);
        jg_plus_i = j_g_plus(i);
        cap_i = cap(i);
        delta_rho = rho_fi - rho_gi;
        
        L_star = Dhi / sqrt(sigma_i / (g * delta_rho));
        N_mu = mu_fi / sqrt(sigma_i * rho_fi * sqrt(sigma_i / (g * delta_rho)));
        r_term = (rho_gi / rho_fi)^(-0.157);
        
        % 初始 alpha
        alpha_old = max(alpha_m(i), 0.0001);
        if alpha_old > 0.99 || isnan(alpha_old)
            alpha_old = 0.0001;
        end
        
        for iter = 1:max_iter
            % === Ozaki 模型的 C0 ===
            C0 = 1.08 - 0.08 * sqrt(rho_gi / rho_fi);
            C0_ozaki_all(i) = C0;
            
            % === Clark vgj 结构复用 ===
            if N_mu <= 0.00225 && L_star <= 30
                vgj_P = 0.0019 * L_star^0.809 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu <= 0.00225 && L_star > 30
                vgj_P = 0.030 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu > 0.00225 && L_star > 30
                vgj_P = 0.92 * r_term * cap_i;
            else
                vgj_P = NaN;
            end
            
            vgj_B = sqrt_term(i) * (1 - alpha_old)^1.75;
            vgj = vgj_B * exp(-1.39 * jg_plus_i) + vgj_P * (1 - exp(-1.39 * jg_plus_i));
            vgj_ozaki_all(i) = vgj;
            
            denom = C0 * jtoti + vgj;
            if denom <= 0
                break;
            end
            
            alpha_new = jgi / denom;
            if ~isfinite(alpha_new) || alpha_new < 0 || alpha_new > 1
                break;
            end
            
            if abs(alpha_new - alpha_old) < tol
                break;
            end
            
            alpha_old = alpha_new;
        end
        
        alpha_ozaki_all(i) = alpha_new;
    end
    
    % 分段误差
    alpha_seg = alpha_ozaki_all(idx_range);
    alpha_m_seg = alpha_m(idx_range);
    valid_idx = alpha_m_seg > 1e-3 & isfinite(alpha_seg);
    rel_err = abs((alpha_seg(valid_idx) - alpha_m_seg(valid_idx)) ./ alpha_m_seg(valid_idx));
    segment_errors_ozaki(k) = mean(rel_err, 'omitnan');
    
    start_idx = start_idx + counts(k);
end

% 加权平均误差
m_rel_OzHi = sum(segment_errors_ozaki .* counts') / total_N;
v_gj_0zHi1 = vgj_ozaki_all;
C0_0zHi1=C0_ozaki_all;


%% Ozaki and Hibiki 模型 (No.2 分段加权)

counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
start_idx = 1;
total_N = sum(counts);

alpha_ozaki_all = NaN(total_N,1);
vgj_ozaki_all   = NaN(total_N,1);
C0_ozaki_all    = NaN(total_N,1);
segment_errors_ozaki = NaN(length(counts),1);

max_iter = 400;
tol = 1e-8;

% 非维量
cap = (sigma .* g .* (rho_f - rho_g) ./ rho_f.^2).^(0.25);
j_g_plus = j_g ./ cap;
sqrt_term = sqrt(2) * cap;

for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;
    for ii = 1:length(idx_range)
        i = idx_range(ii);
        
        % 提取参数
        rho_fi = rho_f(i);
        rho_gi = rho_g(i);
        mu_fi  = mu_f(i);
        jgi = j_g(i);
        jfi = j_f(i);
        jtoti = jgi + jfi;
        Dhi = Dh(i);
        sigma_i = sigma(i);
        jg_plus_i = j_g_plus(i);
        cap_i = cap(i);
        delta_rho = rho_fi - rho_gi;
        
        L_star = Dhi / sqrt(sigma_i / (g * delta_rho));
        N_mu = mu_fi / sqrt(sigma_i * rho_fi * sqrt(sigma_i / (g * delta_rho)));
        r_term = (rho_gi / rho_fi)^(-0.157);
        
        % 初始 alpha
        alpha_old = max(alpha_m(i), 0.0001);
        if alpha_old > 0.99 || isnan(alpha_old)
            alpha_old = 0.0001;
        end
        
        for iter = 1:max_iter
            % === Ozaki 模型的 C0 ===
            C0 = 1.03 - 0.03 * sqrt(rho_gi / rho_fi);
            C0_ozaki_all(i) = C0;
            
            % === Clark vgj 结构复用 ===
            if N_mu <= 0.00225 && L_star <= 30
                vgj_P = 0.0019 * L_star^0.809 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu <= 0.00225 && L_star > 30
                vgj_P = 0.030 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu > 0.00225 && L_star > 30
                vgj_P = 0.92 * r_term * cap_i;
            else
                vgj_P = NaN;
            end
            
            vgj_B = sqrt_term(i) * (1 - alpha_old)^1.75;
            vgj = vgj_B * exp(-1.39 * jg_plus_i) + vgj_P * (1 - exp(-1.39 * jg_plus_i));
            vgj_ozaki_all(i) = vgj;
            
            denom = C0 * jtoti + vgj;
            if denom <= 0
                break;
            end
            
            alpha_new = jgi / denom;
            if ~isfinite(alpha_new) || alpha_new < 0 || alpha_new > 1
                break;
            end
            
            if abs(alpha_new - alpha_old) < tol
                break;
            end
            
            alpha_old = alpha_new;
        end
        
        alpha_ozaki_all(i) = alpha_new;
    end
    % 分段误差
    alpha_seg = alpha_ozaki_all(idx_range);
    alpha_m_seg = alpha_m(idx_range);
    valid_idx = alpha_m_seg > 1e-3 & isfinite(alpha_seg);
    rel_err = abs((alpha_seg(valid_idx) - alpha_m_seg(valid_idx)) ./ alpha_m_seg(valid_idx));
    segment_errors_ozaki(k) = mean(rel_err, 'omitnan');
    start_idx = start_idx + counts(k);
end

% 加权平均误差
m_rel_OzHi2 = sum(segment_errors_ozaki .* counts') / total_N;
v_gj_0zHi2 = vgj_ozaki_all;
C0_0zHi2=C0_ozaki_all;


%% Ren et al. 模型（完整、严谨版本：含 B_sf）
counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
start_idx = 1;
total_N = sum(counts);

alpha_ren_all = NaN(total_N,1);
C0_ren_all    = NaN(total_N,1);
vgj_ren_all   = NaN(total_N,1);
segment_errors_ren = NaN(length(counts),1);

max_iter = 400;
tol = 1e-8;

cap = (sigma .* g .* (rho_f - rho_g) ./ rho_f.^2).^(0.25);

for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;

    for ii = 1:length(idx_range)
        i = idx_range(ii);

        % 提取物理参数
        rho_fi = rho_f(i); rho_gi = rho_g(i);
        jfi = j_f(i); jgi = j_g(i); jtoti = jfi + jgi;
        sigma_i = sigma(i); p_i = p(i); d_i = d(i);
        cap_i = cap(i);

        % 计算 L_max 和 D_b
        Lmax = sqrt(2) * p_i - d_i;
        Db = sqrt(0.4 * sigma_i / (g * (rho_fi - rho_gi)));

        ratio = Db / Lmax;
        if ratio < 0.6
            B_sf = 1 - ratio / 0.9;
        else
            B_sf = 0.120 * (ratio)^(-2);
        end

        % 初始 alpha
        alpha_mi = alpha_m(i);
        alpha_old = max(alpha_mi, 0.01);
        if alpha_old > 0.99 || isnan(alpha_old)
            alpha_old = 0.01;
        end

        for iter = 1:max_iter
            % === C0(alpha) 显式表达 ===
            C0 = 1.14 * alpha_old^(-0.57) * (0.88 - exp(-31.3 * alpha_old));
            C0_ren_all(i) = C0;

            % === vgj(alpha) ===
            vgj = B_sf * sqrt(2) * cap_i * (1 - alpha_old)^1.75;
            vgj_ren_all(i) = vgj;

            % 更新 alpha
            denom = C0 * jtoti + vgj;
            if denom <= 0
                break;
            end

            alpha_new = jgi / denom;
            if ~isfinite(alpha_new) || alpha_new < 0 || alpha_new > 1
                break;
            end

            if abs(alpha_new - alpha_old) < tol
                break;
            end

            alpha_old = alpha_new;
        end

        alpha_ren_all(i) = alpha_new;
    end

    % 分段误差
    alpha_seg = alpha_ren_all(idx_range);
    alpha_m_seg = alpha_m(idx_range);
    valid_idx = alpha_m_seg > 1e-3 & isfinite(alpha_seg);
    rel_err = abs((alpha_seg(valid_idx) - alpha_m_seg(valid_idx)) ./ alpha_m_seg(valid_idx));
    segment_errors_ren(k) = mean(rel_err, 'omitnan');

    start_idx = start_idx + counts(k);
end

% 加权平均误差
m_rel_Ren = sum(segment_errors_ren .* counts') / total_N;
v_gj_Ren=vgj_ren_all;

%% Ye et al. 模型（分段 + 迭代求解 alpha）
counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
start_idx = 1;
total_N = sum(counts);

alpha_ye_all = NaN(total_N,1);
C0_ye_all    = NaN(total_N,1);
vgj_ye_all   = NaN(total_N,1);
segment_errors_ye = NaN(length(counts),1);

max_iter = 400;
tol = 1e-8;

% 非维变量
cap = (sigma .* g .* (rho_f - rho_g) ./ rho_f.^2).^(0.25);
j_g_plus = j_g ./ cap;
j_f_plus = j_f ./ cap;
sqrt_term = sqrt(2) * cap;

for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;
    for ii = 1:length(idx_range)
        i = idx_range(ii);
        
        % 物性与变量
        rho_fi = rho_f(i); rho_gi = rho_g(i);
        mu_fi = mu_f(i); Dhi = Dh(i);
        jgi = j_g(i); jfi = j_f(i); jtoti = jgi + jfi;
        sigma_i = sigma(i);
        cap_i = cap(i);
        jg_plus_i = j_g_plus(i);
        jf_plus_i = j_f_plus(i);
        nu_fi = mu_fi / rho_fi;
        delta_rho = rho_fi - rho_gi;
        
        % 非维数
        L_star = Dhi / sqrt(sigma_i / (g * delta_rho));
        N_mu   = mu_fi / sqrt(sigma_i * rho_fi * sqrt(sigma_i / (g * delta_rho)));
        r_term = (rho_gi / rho_fi)^(-0.157);
        
        % 初始 alpha
        alpha_old = max(alpha_m(i), 0.0001);
        if alpha_old > 0.99 || isnan(alpha_old)
            alpha_old = 0.0001;
        end
        
        for iter = 1:max_iter
            % vgj_P（Clark结构）
            if N_mu <= 0.00225 && L_star <= 30
                vgj_P = 0.0019 * L_star^0.809 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu <= 0.00225 && L_star > 30
                vgj_P = 0.030 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu > 0.00225 && L_star > 30
                vgj_P = 0.92 * r_term * cap_i;
            else
                vgj_P = NaN;
            end
            
            vgj_B = sqrt_term(i) * (1 - alpha_old)^1.75;
            vgj = vgj_B * exp(-1.39 * jg_plus_i) + vgj_P * (1 - exp(-1.39 * jg_plus_i));
            vgj_ye_all(i) = vgj;
            
            % jf_crit 修正版
            jf_crit = 1.92 * Dhi^0.429 * (sigma_i / rho_fi)^0.089 * nu_fi^(-0.072) * (g * delta_rho / rho_fi)^0.446;
            
            % alpha_crit
            if jfi < jf_crit
                alpha_crit = 0.0463 * jf_plus_i + 0.253;
            else
                alpha_crit = min(0.0284 * jf_plus_i + 0.125, 0.52);
            end
            
            % 用临时 C0 估算 jplus_coo_max
            C0_temp = 1.2;
            m = 1 / (1 - alpha_crit * C0_temp);
            b = vgj / cap_i * alpha_crit / (1 - alpha_crit * C0_temp);
            jplus_coo_max = m * jf_plus_i + b;
            
            % 正确计算 C_inf_H(jplus_coo_max)
            C_inf_H_coomax = 1 + 3.45 * exp(-0.52 * jplus_coo_max^0.51);
            C_inf_H_current = 1 + 3.45 * exp(-0.52 * jg_plus_i^0.51);
            
            % 分段计算 C_inf
            if jg_plus_i <= jplus_coo_max
                numer = C_inf_H_coomax - 1;
                denom = jplus_coo_max - jf_plus_i;
                C_inf = (numer / denom) * jg_plus_i + 1;
            else
                C_inf = C_inf_H_current;
            end
            
            % C0 更新
            C0 = C_inf - (C_inf - 1) * sqrt(rho_gi / rho_fi);
            C0_ye_all(i) = C0;
            
            % alpha 更新
            denom = C0 * jtoti + vgj;
            if denom <= 0
                break;
            end
            
            alpha_new = jgi / denom;
            if ~isfinite(alpha_new) || alpha_new < 0 || alpha_new > 1
                break;
            end
            
            if abs(alpha_new - alpha_old) < tol
                break;
            end
            
            alpha_old = alpha_new;
        end
        
        alpha_ye_all(i) = alpha_new;
    end
    
    % 分段误差统计
    alpha_seg = alpha_ye_all(idx_range);
    alpha_m_seg = alpha_m(idx_range);
    valid_idx = alpha_m_seg > 1e-3 & isfinite(alpha_seg);
    rel_err = abs((alpha_seg(valid_idx) - alpha_m_seg(valid_idx)) ./ alpha_m_seg(valid_idx));
    segment_errors_ye(k) = mean(rel_err, 'omitnan');
    
    start_idx = start_idx + counts(k);
end

% 加权平均误差
m_rel_Ye = sum(segment_errors_ye .* counts') / total_N;
v_gj_Ye = vgj_ye_all;

%% Schlegel and Hibiki 模型（分段加权）
counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
start_idx = 1;
total_N = sum(counts);

alpha_schlegel_all = NaN(total_N,1);
C0_schlegel_all    = NaN(total_N,1);
vgj_schlegel_all   = NaN(total_N,1);
segment_errors_schlegel = NaN(length(counts),1);

max_iter = 400;
tol = 1e-8;

% 非维变量
cap = (sigma .* g .* (rho_f - rho_g) ./ rho_f.^2).^(0.25);
j_g_plus = j_g ./ cap;
j_f_plus = j_f ./ cap;
sqrt_term = sqrt(2) * cap;

for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;
    for ii = 1:length(idx_range)
        i = idx_range(ii);
        
        % 物性参数
        rho_fi = rho_f(i); rho_gi = rho_g(i);
        mu_fi = mu_f(i); Dhi = Dh(i);
        jgi = j_g(i); jfi = j_f(i); jtoti = jgi + jfi;
        sigma_i = sigma(i);
        cap_i = cap(i);
        jg_plus_i = j_g_plus(i);
        jf_plus_i = j_f_plus(i);
        delta_rho = rho_fi - rho_gi;
        nu_fi = mu_fi / rho_fi;
        
        % 非维数
        L_star = Dhi / sqrt(sigma_i / (g * delta_rho));
        N_mu   = mu_fi / sqrt(sigma_i * rho_fi * sqrt(sigma_i / (g * delta_rho)));
        r_term = (rho_gi / rho_fi)^(-0.157);
        
        % 初始 alpha
        alpha_old = max(alpha_m(i), 0.0001);
        if alpha_old > 0.99 || isnan(alpha_old)
            alpha_old = 0.0001;
        end
        
        for iter = 1:max_iter
            % === vgj Clark结构 ===
            if N_mu <= 0.00225 && L_star <= 30
                vgj_P = 0.0019 * L_star^0.809 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu <= 0.00225 && L_star > 30
                vgj_P = 0.030 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu > 0.00225 && L_star > 30
                vgj_P = 0.92 * r_term * cap_i;
            else
                vgj_P = NaN;
            end
            
            vgj_B = sqrt_term(i) * (1 - alpha_old)^1.75;
            vgj = vgj_B * exp(-1.39 * jg_plus_i) + vgj_P * (1 - exp(-1.39 * jg_plus_i));
            vgj_schlegel_all(i) = vgj;
            
            % alpha_crit
            alpha_crit = min(0.0284 * jf_plus_i + 0.125, 0.52);
            
            % 临时 C0 用于 jplus_coo_max
            C0_temp = 1.2;
            m = 1 / (1 - alpha_crit * C0_temp);
            b = (vgj / cap_i) * alpha_crit / (1 - alpha_crit * C0_temp);
            jplus_coo_max = m * jf_plus_i + b;
            
            % === F 因子 ===
            F = min(1, max(0, 1.70 - 582 * (rho_gi / rho_fi)));
            
            % C_inf_H 需要用 coo_max 来计算（图中明确指出）
            C_inf_H_coomax = 1.1 + 1.84 * exp(-0.100 * jplus_coo_max) * F;
            C_inf_H_current = 1.1 + 1.84 * exp(-0.100 * jg_plus_i) * F;
            
            % 分段选取 C_inf
            if jg_plus_i <= jplus_coo_max
                numer = C_inf_H_coomax - 1;
                denom = jplus_coo_max - jf_plus_i;
                C_inf = (numer / denom) * jg_plus_i + 1;
            else
                C_inf = C_inf_H_current;
            end
            
            % === 正式计算 C0 ===
            C0 = C_inf - (C_inf - 1) * sqrt(rho_gi / rho_fi);
            C0_schlegel_all(i) = C0;
            
            % alpha 更新
            denom = C0 * jtoti + vgj;
            if denom <= 0
                break;
            end
            
            alpha_new = jgi / denom;
            if ~isfinite(alpha_new) || alpha_new < 0 || alpha_new > 1
                break;
            end
            
            if abs(alpha_new - alpha_old) < tol
                break;
            end
            
            alpha_old = alpha_new;
        end
        
        alpha_schlegel_all(i) = alpha_new;
    end
    
    % 分段误差
    alpha_seg = alpha_schlegel_all(idx_range);
    alpha_m_seg = alpha_m(idx_range);
    valid_idx = alpha_m_seg > 1e-3 & isfinite(alpha_seg);
    rel_err = abs((alpha_seg(valid_idx) - alpha_m_seg(valid_idx)) ./ alpha_m_seg(valid_idx));
    segment_errors_schlegel(k) = mean(rel_err, 'omitnan');
    
    start_idx = start_idx + counts(k);
end

% 总加权平均误差
m_rel_schlegel = sum(segment_errors_schlegel .* counts') / total_N;
m_rel_SH = m_rel_schlegel;




%% Gui et al. 模型

counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
start_idx = 1;
total_N = sum(counts);

alpha_gui_all = NaN(total_N,1);
C0_gui_all    = NaN(total_N,1);
vgj_gui_all   = NaN(total_N,1);
segment_errors_gui = NaN(length(counts),1);

max_iter = 400;
tol = 1e-8;
% 非维变量
cap = (sigma .* g .* (rho_f - rho_g) ./ rho_f.^2).^(0.25);
for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;
    
    for ii = 1:length(idx_range)
        i = idx_range(ii);
        
        % 提取参数
        rho_fi = rho_f(i); rho_gi = rho_g(i);
        jfi = j_f(i); jgi = j_g(i); jtoti = jfi + jgi;
        Dhi = Dh(i); Gi = G(i);
        alpha_mi = alpha_m(i);
        cap_i = cap(i);
        
        
        % Froude 数
        Fr_f = Gi.^2 ./ ( rho_fi.^2 .* g .* Dhi);
        
        % C_inf
        C_inf = 0.80 + 0.34 * exp(-0.36 * Fr_f^(-1.5));
        
        % C0 计算
        C0 = C_inf - (C_inf - 1) * sqrt(rho_gi / rho_fi);
        C0_gui_all(i) = C0;
        
        % 初始 alpha
        alpha_old = max(alpha_mi, 0.01);
        if alpha_old > 0.99 || isnan(alpha_old)
            alpha_old = 0.01;
        end
        
        for iter = 1:max_iter
            % === vgj 条件判断 ===
            if Gi >= 250 && alpha_old > 0.7
                % 修正后的上半段公式
                numerator = 0.5 * (1 - alpha_old);
                denominator = alpha_old - 0.5 * (1 - alpha_old);
                vgj = (numerator / denominator) * ...
                    (67.61 - 33.92 * exp(8.30 * (rho_gi / rho_fi))).* cap_i;
            else
                % 其余情况
                vgj = (4.02 - 330.92 * exp(-1.31 * sqrt(rho_fi / rho_gi))).* cap_i;
            end
            
            
            vgj_gui_all(i) = vgj;
            
            % alpha 更新
            denom = C0 * jtoti + vgj;
            if denom <= 0
                break;
            end
            
            alpha_new = jgi / denom;
            if ~isfinite(alpha_new) || alpha_new < 0 || alpha_new > 1
                break;
            end
            
            if abs(alpha_new - alpha_old) < tol
                break;
            end
            
            alpha_old = alpha_new;
        end
        
        alpha_gui_all(i) = alpha_new;
    end
    
    % 分段误差
    alpha_seg = alpha_gui_all(idx_range);
    alpha_m_seg = alpha_m(idx_range);
    valid_idx = alpha_m_seg > 1e-3 & isfinite(alpha_seg);
    rel_err = abs((alpha_seg(valid_idx) - alpha_m_seg(valid_idx)) ./ alpha_m_seg(valid_idx));
    segment_errors_gui(k) = mean(rel_err, 'omitnan');
    
    start_idx = start_idx + counts(k);
end

% 加权平均误差
m_rel_Gui = sum(segment_errors_gui .* counts') / total_N;
v_gj_Gui=vgj_gui_all;




%% Kinoshita and Hibiki 模型（分段 + 迭代结构）
counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
start_idx = 1;
total_N = sum(counts);

alpha_kinoshita_all = NaN(total_N,1);
C0_kinoshita_all    = NaN(total_N,1);
vgj_kinoshita_all   = NaN(total_N,1);
segment_errors_kinoshita = NaN(length(counts),1);

max_iter = 400;
tol = 1e-8;

% 非维变量
cap = (sigma .* g .* (rho_f - rho_g) ./ rho_f.^2).^(0.25);
j_g_plus = j_g ./ cap;
j_f_plus = j_f ./ cap;
sqrt_term = sqrt(2) * cap;

for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;
    for ii = 1:length(idx_range)
        i = idx_range(ii);
        
        % 物性参数
        rho_fi = rho_f(i); rho_gi = rho_g(i);
        mu_fi = mu_f(i); Dhi = Dh(i);
        jgi = j_g(i); jfi = j_f(i); jtoti = jgi + jfi;
        sigma_i = sigma(i);
        cap_i = cap(i);
        jg_plus_i = j_g_plus(i);
        jf_plus_i = j_f_plus(i);
        delta_rho = rho_fi - rho_gi;
        nu_fi = mu_fi / rho_fi;
        
        % 非维数
        L_star = Dhi / sqrt(sigma_i / (g * delta_rho));
        N_mu   = mu_fi / sqrt(sigma_i * rho_fi * sqrt(sigma_i / (g * delta_rho)));
        r_term = (rho_gi / rho_fi)^(-0.157);
        
        % 初始 alpha
        alpha_old = max(alpha_m(i), 0.0001);
        if alpha_old > 0.99 || isnan(alpha_old)
            alpha_old = 0.0001;
        end
        
        for iter = 1:max_iter
            % vgj_P（Clark结构）
            if N_mu <= 0.00225 && L_star <= 30
                vgj_P = 0.0019 * L_star^0.809 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu <= 0.00225 && L_star > 30
                vgj_P = 0.030 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu > 0.00225 && L_star > 30
                vgj_P = 0.92 * r_term * cap_i;
            else
                vgj_P = NaN;
            end
            
            vgj_B = sqrt_term(i) * (1 - alpha_old)^1.75;
            vgj = vgj_B * exp(-1.39 * jg_plus_i) + vgj_P * (1 - exp(-1.39 * jg_plus_i));
            vgj_kinoshita_all(i) = vgj;
            
            % alpha_crit
            alpha_crit = min(0.0284 * jf_plus_i + 0.125, 0.52);
            
            % 临时 C0 用于估算 jplus_coo_max
            C0_temp = 1.2;
            m = 1 / (1 - alpha_crit * C0_temp);
            b = (vgj / cap_i) * alpha_crit / (1 - alpha_crit * C0_temp);
            jplus_coo_max = m * jf_plus_i + b;
            
            % === F 因子 ===
            F = min(1, max(0, 1.70 - 582 * (rho_gi / rho_fi)));
            
            % C_inf_H（必须用 jplus_coo_max）
            C_inf_H_coomax = 1.1 + 1.84 * exp(-0.100 * jplus_coo_max) * F;
            C_inf_H_current = 1.1 + 1.84 * exp(-0.100 * jg_plus_i) * F;
            
            % 分段选取 C_inf
            if jg_plus_i <= jplus_coo_max
                numer = C_inf_H_coomax - 1;
                denom = jplus_coo_max - jf_plus_i;
                C_inf = (numer / denom) * jg_plus_i + 1;
            else
                C_inf = C_inf_H_current;
            end
            
            % C0 正式计算
            C0 = C_inf - (C_inf - 1) * sqrt(rho_gi / rho_fi);
            C0_kinoshita_all(i) = C0;
            
            % alpha 更新
            denom = C0 * jtoti + vgj;
            if denom <= 0
                break;
            end
            
            alpha_new = jgi / denom;
            if ~isfinite(alpha_new) || alpha_new < 0 || alpha_new > 1
                break;
            end
            
            if abs(alpha_new - alpha_old) < tol
                break;
            end
            
            alpha_old = alpha_new;
        end
        
        alpha_kinoshita_all(i) = alpha_new;
    end
    
    % 分段误差
    alpha_seg = alpha_kinoshita_all(idx_range);
    alpha_m_seg = alpha_m(idx_range);
    valid_idx = alpha_m_seg > 1e-3 & isfinite(alpha_seg);
    rel_err = abs((alpha_seg(valid_idx) - alpha_m_seg(valid_idx)) ./ alpha_m_seg(valid_idx));
    segment_errors_kinoshita(k) = mean(rel_err, 'omitnan');
    
    start_idx = start_idx + counts(k);
end

% 总加权平均误差
m_rel_Kinoshita = sum(segment_errors_kinoshita .* counts') / total_N;



%% Hibiki and Tsukamoto 模型（含 Zhang et al. 漂移项）

counts = [65, 71, 20, 27, 35, 126-1,42,  325-9, 208-10, 26,24-4,  16-1,  43];
sum(counts)
65+71+20+27+35+126-1+325-9+208-10+50-4+16-1+42+43
start_idx = 1;
total_N = sum(counts);

alpha_ozaki_all = NaN(total_N,1);
vgj_ozaki_all   = NaN(total_N,1);
C0_ozaki_all    = NaN(total_N,1);
segment_errors_ozaki = NaN(length(counts),1);

max_iter = 400;
tol = 1e-8;

% 非维量
cap = (sigma .* g .* (rho_f - rho_g) ./ rho_f.^2).^(0.25);
j_g_plus = j_g ./ cap;
sqrt_term = sqrt(2) * cap;

for k = 1:length(counts)
    idx_range = start_idx : start_idx + counts(k) - 1;
    for ii = 1:length(idx_range)
        i = idx_range(ii);
        
        % 提取参数
        rho_fi = rho_f(i);
        rho_gi = rho_g(i);
        mu_fi  = mu_f(i);
        jgi = j_g(i);
        jfi = j_f(i);
        jtoti = jgi + jfi;
        Dhi = Dh(i);
        sigma_i = sigma(i);
        jg_plus_i = j_g_plus(i);
        cap_i = cap(i);
        delta_rho = rho_fi - rho_gi;
        
        L_star = Dhi / sqrt(sigma_i / (g * delta_rho));
        N_mu = mu_fi / sqrt(sigma_i * rho_fi * sqrt(sigma_i / (g * delta_rho)));
        r_term = (rho_gi / rho_fi)^(-0.157);
        
        % 初始 alpha
        alpha_old = max(alpha_m(i), 0.0001);
        if alpha_old > 0.99 || isnan(alpha_old)
            alpha_old = 0.0001;
        end
        
        for iter = 1:max_iter
            % === Ozaki 模型的 C0 ===
            C0 = 1.1 - 0.1 * sqrt(rho_gi / rho_fi);
            C0_ozaki_all(i) = C0;
            
            % === Clark vgj 结构复用 ===
            if N_mu <= 0.00225 && L_star <= 30
                vgj_P = 0.0019 * L_star^0.809 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu <= 0.00225 && L_star > 30
                vgj_P = 0.030 * r_term * N_mu^(-0.562) * cap_i;
            elseif N_mu > 0.00225 && L_star > 30
                vgj_P = 0.92 * r_term * cap_i;
            else
                vgj_P = NaN;
            end
            
            vgj_B = sqrt_term(i) * (1 - alpha_old)^1.75;
            vgj = vgj_B * exp(-1.39 * jg_plus_i) + vgj_P * (1 - exp(-1.39 * jg_plus_i));
            vgj_ozaki_all(i) = vgj;
            
            denom = C0 * jtoti + vgj;
            if denom <= 0
                break;
            end
            
            alpha_new = jgi / denom;
            if ~isfinite(alpha_new) || alpha_new < 0 || alpha_new > 1
                break;
            end
            
            if abs(alpha_new - alpha_old) < tol
                break;
            end
            
            alpha_old = alpha_new;
        end
        
        alpha_ozaki_all(i) = alpha_new;
    end
    % 分段误差
    alpha_seg = alpha_ozaki_all(idx_range);
    alpha_m_seg = alpha_m(idx_range);
    valid_idx = alpha_m_seg > 1e-3 & isfinite(alpha_seg);
    rel_err = abs((alpha_seg(valid_idx) - alpha_m_seg(valid_idx)) ./ alpha_m_seg(valid_idx));
    segment_errors_ozaki(k) = mean(rel_err, 'omitnan');
    start_idx = start_idx + counts(k);
end

% 加权平均误差
m_rel_HT = sum(segment_errors_ozaki .* counts') / total_N;
v_gj_HT = vgj_ozaki_all;
C0_HT=C0_ozaki_all;

% 保存 Hibiki and Tsukamoto 的 alpha，供后续 PIE-DF / PIE-DF-E 专家库识别
% 之前这个模型虽然输出了 m_rel_HT，但没有明确保存 alpha_HT_saved，
% 容易在后面的 modelSpec 收集环节因为命名不一致而漏掉。
alpha_HT_saved = alpha_ozaki_all(:);
alpha_HT_saved = alpha_HT_saved(1:total_N);
alpha_HibikiTsukamoto_saved = alpha_HT_saved;

%% Dix (1971) 模型

% 初始化
sqrt_term = (g .* sigma .* (rho_f - rho_g) ./ rho_f.^2).^(1/4);

% 混合速度
j_mix = j_f + j_g;
j_g_ratio = j_g ./ j_mix;
rho_ratio_01 = (rho_g ./ rho_f).^0.1;

% C0
C0_Dix = j_g_ratio .* (1 + (j_g_ratio - 1) .* rho_ratio_01);

% 漂移速度
v_gj_Dix = 2.9 .* sqrt_term;

% 预测 alpha
alpha_Dix = j_g ./ (C0_Dix .* j_mix + v_gj_Dix);

% 相对误差
m_rel_Dix = mean(abs((alpha_Dix - alpha_m) ./ alpha_m), 'omitnan');

%% Sun et al. / Anklam & Miller 模型

P_crit = 22.09e6;  % Pa, 临界压力
sqrt_term = (g .* sigma .* (rho_f - rho_g) ./ rho_f.^2).^(1/4);

% 混合速度比系数
C0_Sun = 1 ./ (0.82 + 0.18 .* (P ./ P_crit));

% 漂移速度
v_gj_Sun = 1.53 .* sqrt_term;

% alpha 预测
alpha_Sun = j_g ./ (C0_Sun .* j_total + v_gj_Sun);

% 平均相对误差
m_rel_Sun = mean(abs((alpha_Sun - alpha_m) ./ alpha_m), 'omitnan');




%% Hori et al. (1993) 模型

% 将压力单位转换为 MPa
P_MPa = P * 1e-6;

% C0
C0_Hori = 6.76e-3 .* P_MPa + 1.026;

% vgj
v_gj_Hori = (5.1e-3 .* G + 6.91e-2) .* (9.42e-2 .* P_MPa.^2 - 1.99 .* P_MPa + 12.6);

% 预测 alpha
alpha_Hori = j_g ./ (C0_Hori .* j_total + v_gj_Hori);

% 相对误差
m_rel_Hori = mean(abs((alpha_Hori - alpha_m) ./ alpha_m), 'omitnan');

%% 显示结果
fprintf('Jowit 平均相对误差 m_rel = %.3f\n', m_rel_Jowitt);
fprintf('Morooka    平均相对误差 m_rel = %.3f\n', m_rel_Morooka);
fprintf('Bestion     平均相对误差 m_rel = %.3f\n', m_rel_Bestion);
fprintf('Chexal-Lellouche 平均相对误差 m_rel = %.3f\n', m_rel_chexal);
fprintf('Maier 平均相对误差 m_rel = %.3f\n', m_rel_MC);
fprintf('Julia et al. 平均相对误差 m_rel = %.3f\n', m_rel_Julia);
fprintf('Paranjape 平均相对误差 m_rel = %.3f\n', m_rel_Paranjape);
fprintf('Kamei 平均相对误差 m_rel = %.3f\n', m_rel_Kamei);
fprintf('Chen 平均相对误差 m_rel = %.3f\n', m_rel_Chen);
fprintf('Ozaki 平均相对误差 m_rel = %.3f\n', m_rel_Ozaki);
fprintf('Clark 平均相对误差 m_rel = %.3f\n', m_rel_Clark);
fprintf('Ozaki and Hibiki (No.1) 平均相对误差 m_rel = %.3f\n', m_rel_OzHi);
fprintf('Ozaki and Hibiki (No.2) 平均相对误差 m_rel = %.3f\n', m_rel_OzHi2);
fprintf('Ren et al. 平均相对误差 m_rel = %.3f\n', m_rel_Ren);
fprintf('Ye et al. 平均相对误差 m_rel = %.3f\n', m_rel_Ye);
fprintf('Schlegel and Hibiki 平均相对误差 m_rel = %.3f\n', m_rel_SH);
fprintf('Gui et al. 平均相对误差 m_rel = %.3f\n', m_rel_Gui);
fprintf('Kinoshita and Hibiki. 平均相对误差 m_rel = %.3f\n', m_rel_Kinoshita);
fprintf('Hibiki and Tsukamoto. 平均相对误差 m_rel = %.3f\n', m_rel_HT);
fprintf('Dix  平均相对误差 m_rel = %.4f\n', m_rel_Dix);
fprintf('Sun / Anklam & Miller 平均相对误差 m_rel = %.3f\n', m_rel_Sun);
fprintf('Hori et al.  平均相对误差 m_rel = %.3f\n', m_rel_Hori);
fprintf('Our model.  平均相对误差 m_rel = %.4f\n', (0.0805*248+0.0924*755)/(248+755));


%% ========================================================================
%  ONE-CLICK APPEND CODE
%  PIE-DF / PIE-DF-E 机器学习嵌入模型 + 13数据库 ML 作图
%
%  使用方法：
%  1. 把本段代码完整粘贴到你原始 1003 点漂移流模型代码最后面；
%  2. 直接运行整个 .m 文件；
%  3. 输出文件夹：ML_Embedded_13DB_Plots
%
%  模型逻辑：
%       PIE-DF   : alpha_PIE  = sum_i w_i(x) * alpha_i
%       PIE-DF-E : alpha_PIEE = alpha_PIE + err(x)
%
%  只画 ML 模型，不画传统模型。
% ========================================================================

fprintf('\n\n============================================================\n');
fprintf('开始运行 PIE-DF / PIE-DF-E 一键机器学习嵌入模型\n');
fprintf('============================================================\n');

rng(2026);

%% ========================================================================
%  0. 基本设置
% ========================================================================

db_counts = [65 71 20 27 35 125 42 316 198 26 20 15 43];
N_total_db = sum(db_counts);

if ~exist('alpha_m','var')
    error('未找到 alpha_m。请确认本代码已经粘贴在原始漂移流模型代码后面。');
end

requiredVars = {'j_f','j_g','rho_f','rho_g','mu_f','mu_g','Re_f','Re_g','d','p','Dh','P','G','sigma'};

for kk = 1:numel(requiredVars)
    if ~exist(requiredVars{kk},'var')
        error('缺少变量 %s。请确认原始模型代码已经正常运行。', requiredVars{kk});
    end
end

if length(alpha_m) < N_total_db
    error('alpha_m 长度为 %d，但13个数据库总点数为 %d。', length(alpha_m), N_total_db);
end

alpha_m = alpha_m(:);
alpha_m = alpha_m(1:N_total_db);

j_f   = j_f(:);     j_f   = j_f(1:N_total_db);
j_g   = j_g(:);     j_g   = j_g(1:N_total_db);
rho_f = rho_f(:);   rho_f = rho_f(1:N_total_db);
rho_g = rho_g(:);   rho_g = rho_g(1:N_total_db);
mu_f  = mu_f(:);    mu_f  = mu_f(1:N_total_db);
mu_g  = mu_g(:);    mu_g  = mu_g(1:N_total_db);
Re_f  = Re_f(:);    Re_f  = Re_f(1:N_total_db);
Re_g  = Re_g(:);    Re_g  = Re_g(1:N_total_db);
d     = d(:);       d     = d(1:N_total_db);
p     = p(:);       p     = p(1:N_total_db);
Dh    = Dh(:);      Dh    = Dh(1:N_total_db);
P     = P(:);       P     = P(1:N_total_db);
G     = G(:);       G     = G(1:N_total_db);
sigma = sigma(:);   sigma = sigma(1:N_total_db);

j_total = j_f + j_g;
g = 9.81;

db_start = zeros(length(db_counts),1);
db_end   = zeros(length(db_counts),1);

s = 1;
for db = 1:length(db_counts)
    db_start(db) = s;
    db_end(db) = s + db_counts(db) - 1;
    s = db_end(db) + 1;
end

fprintf('\n================ 13个数据库分段信息 ================\n');
for db = 1:length(db_counts)
    fprintf('DB %2d: %4d ~ %4d, 点数 = %d\n', ...
        db, db_start(db), db_end(db), db_counts(db));
end

%% ========================================================================
%  1. 补充保存 / 重构部分模型 alpha
% ========================================================================

% Bestion
if exist('alpha_Bestion_all','var') && ~exist('alpha_Bestion_saved','var')
    tmp = alpha_Bestion_all(:);
    if length(tmp) >= N_total_db
        alpha_Bestion_saved = tmp(1:N_total_db);
    end
end

% Chexal
if exist('alpha_chexal_all','var') && ~exist('alpha_Chexal_saved','var')
    tmp = alpha_chexal_all(:);
    if length(tmp) >= N_total_db
        alpha_Chexal_saved = tmp(1:N_total_db);
    end
end

% Chen
if exist('alpha_chen_all','var') && ~exist('alpha_Chen_saved','var')
    tmp = alpha_chen_all(:);
    if length(tmp) >= N_total_db
        alpha_Chen_saved = tmp(1:N_total_db);
    end
end

% Ren
if exist('alpha_ren_all','var') && ~exist('alpha_Ren_saved','var')
    tmp = alpha_ren_all(:);
    if length(tmp) >= N_total_db
        alpha_Ren_saved = tmp(1:N_total_db);
    end
end

% Ye
if exist('alpha_ye_all','var') && ~exist('alpha_Ye_saved','var')
    tmp = alpha_ye_all(:);
    if length(tmp) >= N_total_db
        alpha_Ye_saved = tmp(1:N_total_db);
    end
end

% Schlegel-Hibiki
if exist('alpha_schlegel_all','var') && ~exist('alpha_SH_saved','var')
    tmp = alpha_schlegel_all(:);
    if length(tmp) >= N_total_db
        alpha_SH_saved = tmp(1:N_total_db);
    end
end

% Gui
if exist('alpha_gui_all','var') && ~exist('alpha_Gui_saved','var')
    tmp = alpha_gui_all(:);
    if length(tmp) >= N_total_db
        alpha_Gui_saved = tmp(1:N_total_db);
    end
end

% Kinoshita
if exist('alpha_kinoshita_all','var') && ~exist('alpha_Kinoshita_saved','var')
    tmp = alpha_kinoshita_all(:);
    if length(tmp) >= N_total_db
        alpha_Kinoshita_saved = tmp(1:N_total_db);
    end
end

% Clark：你的 Clark 模型通常存在 alpha_all_model 里
if exist('alpha_all_model','var') && ~exist('alpha_Clark_saved','var')
    tmp = alpha_all_model(:);
    if length(tmp) >= N_total_db
        alpha_Clark_saved = tmp(1:N_total_db);
    end
end

% Ozaki-Hibiki No.1
if exist('C0_0zHi1','var') && exist('v_gj_0zHi1','var') && ~exist('alpha_OzHi1_saved','var')
    C0_tmp = C0_0zHi1(:);
    vj_tmp = v_gj_0zHi1(:);
    if length(C0_tmp) >= N_total_db && length(vj_tmp) >= N_total_db
        alpha_OzHi1_saved = j_g ./ (C0_tmp(1:N_total_db) .* j_total + vj_tmp(1:N_total_db));
    end
end

% Ozaki-Hibiki No.2
if exist('C0_0zHi2','var') && exist('v_gj_0zHi2','var') && ~exist('alpha_OzHi2_saved','var')
    C0_tmp = C0_0zHi2(:);
    vj_tmp = v_gj_0zHi2(:);
    if length(C0_tmp) >= N_total_db && length(vj_tmp) >= N_total_db
        alpha_OzHi2_saved = j_g ./ (C0_tmp(1:N_total_db) .* j_total + vj_tmp(1:N_total_db));
    end
end

% Hibiki and Tsukamoto / Zhang-H-T
% 统一命名，确保该模型能够进入 PIE-DF / PIE-DF-E 专家库。
% 优先使用前面模型段直接迭代得到的 alpha_HT_saved；如果不存在，再用 C0_HT 和 v_gj_HT 重构。
if exist('alpha_HT_saved','var')
    tmp = alpha_HT_saved(:);
    if length(tmp) >= N_total_db
        alpha_HT_saved = tmp(1:N_total_db);
        alpha_HibikiTsukamoto_saved = alpha_HT_saved;
        alpha_Zhang_saved = alpha_HT_saved;
    end
elseif exist('alpha_HibikiTsukamoto_saved','var')
    tmp = alpha_HibikiTsukamoto_saved(:);
    if length(tmp) >= N_total_db
        alpha_HibikiTsukamoto_saved = tmp(1:N_total_db);
        alpha_HT_saved = alpha_HibikiTsukamoto_saved;
        alpha_Zhang_saved = alpha_HibikiTsukamoto_saved;
    end
elseif exist('C0_HT','var') && exist('v_gj_HT','var')
    C0_tmp = C0_HT(:);
    vj_tmp = v_gj_HT(:);
    if length(C0_tmp) >= N_total_db && length(vj_tmp) >= N_total_db
        alpha_HT_saved = j_g ./ (C0_tmp(1:N_total_db) .* j_total + vj_tmp(1:N_total_db));
        alpha_HibikiTsukamoto_saved = alpha_HT_saved;
        alpha_Zhang_saved = alpha_HT_saved;
    end
end

% Hori
if exist('C0_Hori','var') && exist('v_gj_Hori','var') && ~exist('alpha_Hori','var')
    alpha_Hori = j_g ./ (C0_Hori(:) .* j_total + v_gj_Hori(:));
end

% Julia 已在 Julia 模型段结束后立即保存为 alpha_Julia_saved。
% 这里不再从 alpha_all_model 取 Julia，避免被后面的 Clark 结果覆盖。

%% ========================================================================
%  2. 收集传统模型 alpha 矩阵 A
% ========================================================================

modelSpec = {
    'Sun et al.',             {'alpha_Sun_saved','alpha_Sun','alpha_sun_all'}
    'Jowitt',                 {'alpha_Jowitt','alpha_jowitt_all'}
    'Bestion',                {'alpha_Bestion_saved','alpha_Bestion_all','alpha_Bestion'}
    'Morooka',                {'alpha_Morooka','alpha_morooka_all'}
    'Chexal et al.',          {'alpha_Chexal_saved','alpha_chexal_all','alpha_Chexal'}
    'Maier et al.',           {'alpha_MC','alpha_Maier','alpha_maier_all'}
    'Paranjape et al.',       {'alpha_Paranjape','alpha_paranjape_all'}
    'Julia et al.',           {'alpha_Julia_saved','alpha_Julia','alpha_julia_all'}
    'Kamei et al.',           {'alpha_Kamei','alpha_kamei_all'}
    'Chen et al.',            {'alpha_Chen_saved','alpha_chen_all','alpha_Chen'}
    'Clark et al.',           {'alpha_Clark_saved','alpha_Clark','alpha_clark_all'}
    'Ozaki-Hibiki No.1',      {'alpha_OzHi1_saved'}
    'Ozaki-Hibiki No.2',      {'alpha_OzHi2_saved'}
    'Ren et al.',             {'alpha_Ren_saved','alpha_ren_all','alpha_Ren'}
    'Ye et al.',              {'alpha_Ye_saved','alpha_ye_all','alpha_Ye'}
    'Schlegel-Hibiki',        {'alpha_SH_saved','alpha_schlegel_all','alpha_SH'}
    'Gui et al.',             {'alpha_Gui_saved','alpha_gui_all','alpha_Gui'}
    'Kinoshita et al.',       {'alpha_Kinoshita_saved','alpha_kinoshita_all','alpha_Kinoshita'}
    'Hibiki-Tsukamoto',       {'alpha_HT_saved','alpha_HibikiTsukamoto_saved','alpha_Zhang_saved','alpha_HT','alpha_zhang_all','alpha_Zhang'}
    'Dix',                    {'alpha_Dix'}
    'Hori et al.',            {'alpha_Hori'}
    };

% ------------------------------------------------------------------------
%  不参与 PIE-DF / PIE-DF-E 的传统模型
%  说明：这里只从机器学习嵌入模型的专家库 A 中剔除，
%  前面各模型自身的 alpha、m_rel 计算和输出不删除。
% ------------------------------------------------------------------------
excludeFromPIEDF = { ...
    'Hori et al.', ...
    'Schlegel-Hibiki', ...
    'Ozaki-Hibiki No.1', ...
    'Ozaki', ...
    'Dix'};

keepModelSpec = true(size(modelSpec,1),1);
for ee = 1:numel(excludeFromPIEDF)
    keepModelSpec = keepModelSpec & ~strcmpi(modelSpec(:,1), excludeFromPIEDF{ee});
end

excludedActuallyFound = modelSpec(~keepModelSpec,1);
modelSpec = modelSpec(keepModelSpec,:);

fprintf('\n================ PIE-DF / PIE-DF-E 专家库剔除模型 ================\n');
if isempty(excludedActuallyFound)
    fprintf('未在 modelSpec 中找到需要剔除的模型名称。\n');
else
    for ee = 1:numel(excludedActuallyFound)
        fprintf('已剔除，不参与 PIE-DF 和 PIE-DF-E：%s\n', excludedActuallyFound{ee});
    end
end
fprintf('保留并参与 PIE-DF / PIE-DF-E：Hibiki-Tsukamoto，候选变量 alpha_HT_saved / alpha_HibikiTsukamoto_saved / alpha_Zhang_saved\n');

A = [];
modelNames = {};
usedVarNames = {};
missingModels = {};
shortModels = {};

fprintf('\n================ 模型 alpha 向量收集结果 ================\n');

for ii = 1:size(modelSpec,1)

    showName = modelSpec{ii,1};
    candidateVars = modelSpec{ii,2};

    found = false;
    foundShort = false;
    shortInfo = '';

    for jj = 1:numel(candidateVars)

        vname = candidateVars{jj};

        if exist(vname,'var')

            tmp = eval(vname);

            if isnumeric(tmp)

                tmp = tmp(:);

                if length(tmp) >= N_total_db

                    tmp = tmp(1:N_total_db);

                    validRatio = mean(isfinite(tmp) & tmp > 0 & tmp < 1);

                    if validRatio >= 0.60

                        duplicate = false;

                        for mm = 1:size(A,2)
                            diffCheck = A(:,mm) - tmp;
                            sameRatio = mean(abs(diffCheck) < 1e-12 | ~isfinite(diffCheck));
                            if sameRatio > 0.999
                                duplicate = true;
                                break;
                            end
                        end

                        if ~duplicate
                            A = [A, tmp];
                            modelNames{end+1,1} = showName;
                            usedVarNames{end+1,1} = vname;

                            fprintf('已收集：%-22s  使用变量：%-28s  有效率：%.1f%%\n', ...
                                showName, vname, 100*validRatio);
                        else
                            fprintf('跳过重复：%-22s  变量：%-28s\n', showName, vname);
                        end

                        found = true;
                        break;

                    else
                        fprintf('跳过：%-22s  变量：%-28s  有效率过低：%.1f%%\n', ...
                            showName, vname, 100*validRatio);
                    end

                else
                    foundShort = true;
                    shortInfo = sprintf('%s 长度=%d', vname, length(tmp));
                end
            end
        end
    end

    if ~found
        missingModels{end+1,1} = showName;
        if foundShort
            shortModels{end+1,1} = sprintf('%s: %s', showName, shortInfo);
            fprintf('未收集：%-22s  原因：长度不足1003；%s\n', showName, shortInfo);
        else
            fprintf('未收集：%-22s  原因：未找到对应 alpha 变量\n', showName);
        end
    end
end

M = size(A,2);

fprintf('\n成功进入 ML 的传统模型数 M = %d\n', M);

if M < 2
    error('可用传统模型少于2个，无法训练 PIE-DF。');
end

if ~isempty(missingModels)
    fprintf('\n未进入 ML 的模型：\n');
    for ii = 1:numel(missingModels)
        fprintf('  - %s\n', missingModels{ii});
    end
end

%% ========================================================================
%  3. 清洗传统模型预测矩阵
% ========================================================================

for mm = 1:M

    col = A(:,mm);

    good = isfinite(col) & col > 0 & col < 1;

    medVal = median(col(good), 'omitnan');

    if ~isfinite(medVal)
        medVal = 0.5;
    end

    col(~good) = medVal;
    col = min(max(col, 1e-5), 0.999);

    A(:,mm) = col;
end

Y = alpha_m(:);
Y = min(max(Y, 1e-5), 0.999);

%% ========================================================================
%  4. 构造 ML 输入特征
% ========================================================================

delta_rho = max(rho_f - rho_g, eps);
cap_ml = (sigma .* g .* delta_rho ./ rho_f.^2).^(1/4);

rhoRatio = rho_g ./ rho_f;
muRatio  = mu_g ./ mu_f;

pitchRatio = p ./ d;
gapRatio = (p - d) ./ d;
areaLike = p.^2 - pi .* d.^2 ./ 4;

Fr_f_ml = j_f ./ sqrt(g .* Dh);
Fr_g_ml = j_g ./ sqrt(g .* Dh);

We_f = rho_f .* j_f.^2 .* Dh ./ sigma;
We_g = rho_g .* j_g.^2 .* Dh ./ sigma;

La = sqrt(sigma ./ (g .* delta_rho));
Lstar = Dh ./ La;

Nmu = mu_f ./ sqrt(sigma .* rho_f .* La);

jgPlus = j_g ./ cap_ml;
jfPlus = j_f ./ cap_ml;
jPlus  = j_total ./ cap_ml;

P_MPa = P * 1e-6;

Xcont = [
    j_f, j_g, j_total, ...
    P_MPa, G, ...
    rho_f, rho_g, rhoRatio, ...
    mu_f, mu_g, muRatio, ...
    Re_f, Re_g, ...
    d, p, Dh, pitchRatio, gapRatio, areaLike, ...
    sigma, cap_ml, ...
    Fr_f_ml, Fr_g_ml, We_f, We_g, ...
    Lstar, Nmu, ...
    jfPlus, jgPlus, jPlus
    ];

featureNames = {
    'j_f','j_g','j_total', ...
    'P_MPa','G', ...
    'rho_f','rho_g','rho_g_over_rho_f', ...
    'mu_f','mu_g','mu_g_over_mu_f', ...
    'Re_f','Re_g', ...
    'd','p','Dh','p_over_d','gap_over_d','area_like', ...
    'sigma','cap', ...
    'Fr_f','Fr_g','We_f','We_g', ...
    'Lstar','Nmu', ...
    'jf_plus','jg_plus','j_plus'
    };

groupID = zeros(N_total_db,1);

for db = 1:length(db_counts)
    idx = db_start(db):db_end(db);
    groupID(idx) = db;
end

numGroups = max(groupID);

Xgroup = zeros(N_total_db, numGroups);

for ii = 1:N_total_db
    Xgroup(ii, groupID(ii)) = 1;
end

validRows = isfinite(Y) & Y > 1e-5 & Y < 1 & ...
            all(isfinite(Xcont),2) & ...
            all(isfinite(A),2);

fprintf('\n有效样本数量 = %d / %d\n', sum(validRows), N_total_db);

Xcont_valid  = Xcont(validRows,:);
Xgroup_valid = Xgroup(validRows,:);
Y_valid      = Y(validRows);
A_valid      = A(validRows,:);
groupID_valid = groupID(validRows);

Nvalid = length(Y_valid);

%% ========================================================================
%  5. 构造监督权重 Wtrue
% ========================================================================

absErr = abs(A_valid - Y_valid);

tau = median(absErr(:), 'omitnan');

if ~isfinite(tau) || tau <= 0
    tau = 0.05;
end

tau = max(tau, 0.02);

Wtrue = exp(-absErr ./ tau);
Wtrue = Wtrue ./ sum(Wtrue, 2);

[~, bestModelEachPoint] = max(Wtrue, [], 2);

fprintf('\n================ 逐点最优传统模型频率 ================\n');

for mm = 1:M
    fprintf('%-22s : %.2f %%\n', modelNames{mm}, 100*mean(bestModelEachPoint == mm));
end

%% ========================================================================
%  6. 按数据库分层划分训练 / 验证 / 测试
% ========================================================================

idxAll = (1:Nvalid)';

idxTrain = [];
idxVal   = [];
idxTest  = [];

for db = 1:max(groupID_valid)

    idxK = idxAll(groupID_valid == db);
    idxK = idxK(randperm(numel(idxK)));

    nK = numel(idxK);

    if nK <= 2
        idxTrain = [idxTrain; idxK];
        continue;
    end

    nTrain = max(1, round(0.70*nK));
    nVal   = max(1, round(0.15*nK));
    nTest  = nK - nTrain - nVal;

    if nTest < 1
        nTrain = max(1, nTrain - 1);
        nTest = 1;
    end

    idxTrain = [idxTrain; idxK(1:nTrain)];
    idxVal   = [idxVal; idxK(nTrain+1:nTrain+nVal)];
    idxTest  = [idxTest; idxK(nTrain+nVal+1:end)];
end

idxTrain = idxTrain(randperm(numel(idxTrain)));
idxVal   = idxVal(randperm(numel(idxVal)));
idxTest  = idxTest(randperm(numel(idxTest)));

if isempty(idxVal)
    idxVal = idxTest;
end

fprintf('\n训练集: %d, 验证集: %d, 测试集: %d\n', ...
    numel(idxTrain), numel(idxVal), numel(idxTest));

%% ========================================================================
%  7. 标准化输入
% ========================================================================

muX = mean(Xcont_valid(idxTrain,:), 1, 'omitnan');
stdX = std(Xcont_valid(idxTrain,:), 0, 1, 'omitnan');

stdX(stdX < 1e-12 | ~isfinite(stdX)) = 1;

XcontN_valid = (Xcont_valid - muX) ./ stdX;
X_valid = [XcontN_valid, Xgroup_valid];

numFeatures = size(X_valid,2);

%% ========================================================================
%  8. 单个传统模型测试集误差
% ========================================================================

fprintf('\n================ 单个传统模型测试集误差 ================\n');

singleMetrics = table( ...
    strings(M,1), ...
    zeros(M,1), ...
    zeros(M,1), ...
    zeros(M,1), ...
    zeros(M,1), ...
    'VariableNames', {'Model','ARE','MAE','RMSE','R2'} );

for mm = 1:M

    met = calcMetricsPIEOneClick_R2FIX20260513(Y_valid(idxTest), A_valid(idxTest,mm));

    fprintf('%-22s ARE=%.4f  MAE=%.4f  RMSE=%.4f  R2=%.4f\n', ...
        modelNames{mm}, met.ARE, met.MAE, met.RMSE, met.R2);

    singleMetrics.Model(mm) = string(modelNames{mm});
    singleMetrics.ARE(mm)   = met.ARE;
    singleMetrics.MAE(mm)   = met.MAE;
    singleMetrics.RMSE(mm)  = met.RMSE;
    singleMetrics.R2(mm)    = met.R2;
end

%% ========================================================================
%  9-10. 多随机种子重复训练：PIE-DF + PIE-DF-E
%
%  审稿人关注点：单次训练不能证明对训练随机性的稳健性。
%  本版本采用“固定数据划分 + 改变训练随机种子”的设计：
%    1) 第6节仍使用最前面的 rng(2026) 产生固定 train/val/test split；
%    2) 下面仅改变训练随机种子，因此每次实验使用完全相同的输入、
%       完全相同的训练/验证/测试样本；
%    3) 随机种子会影响网络初始化、dropout 和每个 epoch 的 mini-batch shuffle；
%    4) 自动输出每个 seed 的 Test ARE/MAE/RMSE/R2，以及 mean/std/min/max；
%    5) seed = 2026 的网络被保存并恢复，供后续第11节开始的全部作图、
%       1003点预测和权重分析继续使用，从而保持原稿主结果口径不变。
%
%  默认使用10个训练随机种子：2026~2035。
%  如需更多重复，只改 trainingSeeds 即可，例如：
%       trainingSeeds = 2026:2045;   % 20次
% ========================================================================

trainingSeeds = 2026:2035;
representativeSeed = 2026;
nSeeds = numel(trainingSeeds);

fprintf('\n\n============================================================\n');
fprintf('开始多随机种子稳健性实验：固定 split，仅改变训练随机性\n');
fprintf('固定数据划分 seed = 2026；训练 seeds = %d ~ %d；共 %d 次\n', ...
    trainingSeeds(1), trainingSeeds(end), nSeeds);
fprintf('============================================================\n');

% ------------------------------------------------------------------------
% 固定网络结构与超参数
% ------------------------------------------------------------------------
layersExpertTemplate = [
    featureInputLayer(numFeatures, "Normalization","none", "Name","input")
    fullyConnectedLayer(96, "Name","fc1")
    tanhLayer("Name","tanh1")
    dropoutLayer(0.05, "Name","drop1")
    fullyConnectedLayer(64, "Name","fc2")
    tanhLayer("Name","tanh2")
    fullyConnectedLayer(32, "Name","fc3")
    tanhLayer("Name","tanh3")
    fullyConnectedLayer(M, "Name","logits")
    ];

layersErrTemplate = [
    featureInputLayer(numFeatures, "Normalization","none", "Name","input")
    fullyConnectedLayer(48, "Name","fc1")
    tanhLayer("Name","tanh1")
    dropoutLayer(0.05, "Name","drop1")
    fullyConnectedLayer(24, "Name","fc2")
    tanhLayer("Name","tanh2")
    fullyConnectedLayer(1, "Name","err")
    ];

numEpochs = 800;
numEpochsErr = 400;

miniBatchSize = min(64, max(16, floor(numel(idxTrain)/3)));
miniBatchSizeErr = miniBatchSize;

learnRate0 = 8e-4;
learnRateErr0 = 5e-4;

lambdaCE = 0.30;
lambdaL2W = 1e-4;
lambdaErrL2 = 1e-4;

gradDecay = 0.9;
sqGradDecay = 0.999;

% ------------------------------------------------------------------------
% 用同一个固定测试集寻找“ARE最优传统模型”，用于判断提升是否超过随机波动
% ------------------------------------------------------------------------
validBaseARE = isfinite(singleMetrics.ARE);
if any(validBaseARE)
    baseRows = find(validBaseARE);
    [bestTraditionalARE, iiBase] = min(singleMetrics.ARE(validBaseARE));
    bestTraditionalIdx = baseRows(iiBase);
    bestTraditionalModelARE = string(singleMetrics.Model(bestTraditionalIdx));
else
    bestTraditionalARE = NaN;
    bestTraditionalIdx = NaN;
    bestTraditionalModelARE = "N/A";
end

fprintf('\n固定测试集上 ARE 最优的单一传统模型：%s，ARE = %.6f\n', ...
    char(bestTraditionalModelARE), bestTraditionalARE);

% ------------------------------------------------------------------------
% 逐 seed 结果表
% ------------------------------------------------------------------------
seedResults = table( ...
    trainingSeeds(:), ...
    NaN(nSeeds,1), NaN(nSeeds,1), NaN(nSeeds,1), NaN(nSeeds,1), ...
    NaN(nSeeds,1), NaN(nSeeds,1), NaN(nSeeds,1), NaN(nSeeds,1), ...
    NaN(nSeeds,1), NaN(nSeeds,1), ...
    'VariableNames', { ...
    'Seed', ...
    'PIE_DF_Test_ARE','PIE_DF_Test_MAE','PIE_DF_Test_RMSE','PIE_DF_Test_R2', ...
    'PIE_DF_E_Test_ARE','PIE_DF_E_Test_MAE','PIE_DF_E_Test_RMSE','PIE_DF_E_Test_R2', ...
    'PIE_DF_ARE_Improvement_vs_BestTraditional_percent', ...
    'PIE_DF_E_ARE_Improvement_vs_BestTraditional_percent'} );

haveRepresentative = false;

% 代表性 seed=2026 的网络/历史/指标，用于后续原有流程
netExpertRepresentative = [];
netErrRepresentative = [];
historyExpertRepresentative = [];
historyErrRepresentative = [];

fprintf('\n训练集=%d，验证集=%d，测试集=%d（所有 seeds 完全一致）\n', ...
    numel(idxTrain), numel(idxVal), numel(idxTest));

% ========================================================================
% 多随机种子主循环
% ========================================================================
for iSeed = 1:nSeeds

    thisSeed = trainingSeeds(iSeed);

    % 关键：只在网络训练开始前重置随机种子。
    % 固定的 idxTrain / idxVal / idxTest 不再重新划分。
    rng(thisSeed, 'twister');

    fprintf('\n\n############################################################\n');
    fprintf('Multi-seed run %d/%d | training seed = %d\n', iSeed, nSeeds, thisSeed);
    fprintf('############################################################\n');

    %% --------------------------------------------------------------------
    %  9A. 训练 PIE-DF 权重网络
    % ---------------------------------------------------------------------
    layersExpert = layersExpertTemplate;
    netExpert = dlnetwork(layerGraph(layersExpert));

    learnRate = learnRate0;
    avgGrad = [];
    avgSqGrad = [];
    iteration = 0;

    historyExpert = zeros(numEpochs, 4);

    for epoch = 1:numEpochs

        % 受 thisSeed 控制的 mini-batch 顺序
        idxTrainShuffle = idxTrain(randperm(numel(idxTrain)));

        for start = 1:miniBatchSize:numel(idxTrainShuffle)

            iteration = iteration + 1;

            batchIdx = idxTrainShuffle(start:min(start+miniBatchSize-1, numel(idxTrainShuffle)));

            dlX = dlarray(single(X_valid(batchIdx,:)'), 'CB');
            dlA = dlarray(single(A_valid(batchIdx,:)'), 'CB');
            dlY = dlarray(single(Y_valid(batchIdx)'), 'CB');
            dlWtrue = dlarray(single(Wtrue(batchIdx,:)'), 'CB');

            [loss, gradients] = dlfeval(@expertGradientsPIEOneClick, ...
                netExpert, dlX, dlA, dlY, dlWtrue, lambdaCE, lambdaL2W);

            [netExpert, avgGrad, avgSqGrad] = adamupdate( ...
                netExpert, gradients, avgGrad, avgSqGrad, ...
                iteration, learnRate, gradDecay, sqGradDecay);
        end

        if mod(epoch, 300) == 0
            learnRate = learnRate * 0.5;
        end

        alphaTrainPIE = predictPIEOneClick(netExpert, X_valid(idxTrain,:), A_valid(idxTrain,:));
        alphaValPIE   = predictPIEOneClick(netExpert, X_valid(idxVal,:), A_valid(idxVal,:));

        metTrain = calcMetricsPIEOneClick_R2FIX20260513(Y_valid(idxTrain), alphaTrainPIE);
        metVal   = calcMetricsPIEOneClick_R2FIX20260513(Y_valid(idxVal), alphaValPIE);

        historyExpert(epoch,:) = [double(gather(extractdata(loss))), ...
                                  metVal.MSE, metVal.ARE, metVal.RMSE];

        if mod(epoch, 100) == 0 || epoch == 1 || epoch == numEpochs
            fprintf('[Seed %d][PIE-DF] Epoch %4d/%4d | Train ARE %.4f | Val ARE %.4f | Val RMSE %.4f\n', ...
                thisSeed, epoch, numEpochs, metTrain.ARE, metVal.ARE, metVal.RMSE);
        end
    end

    alphaTrainPIE = predictPIEOneClick(netExpert, X_valid(idxTrain,:), A_valid(idxTrain,:));
    alphaValPIE   = predictPIEOneClick(netExpert, X_valid(idxVal,:), A_valid(idxVal,:));
    alphaTestPIE  = predictPIEOneClick(netExpert, X_valid(idxTest,:), A_valid(idxTest,:));

    metPIETrain = calcMetricsPIEOneClick_R2FIX20260513(Y_valid(idxTrain), alphaTrainPIE);
    metPIEVal   = calcMetricsPIEOneClick_R2FIX20260513(Y_valid(idxVal), alphaValPIE);
    metPIETest  = calcMetricsPIEOneClick_R2FIX20260513(Y_valid(idxTest), alphaTestPIE);

    %% --------------------------------------------------------------------
    %  10A. 训练 PIE-DF-E 残差网络
    % ---------------------------------------------------------------------
    layersErr = layersErrTemplate;
    netErr = dlnetwork(layerGraph(layersErr));

    learnRateErr = learnRateErr0;
    avgGrad = [];
    avgSqGrad = [];
    iteration = 0;

    historyErr = zeros(numEpochsErr, 3);

    for epoch = 1:numEpochsErr

        % 受同一个 thisSeed 的随机流控制
        idxTrainShuffle = idxTrain(randperm(numel(idxTrain)));

        for start = 1:miniBatchSizeErr:numel(idxTrainShuffle)

            iteration = iteration + 1;

            batchIdx = idxTrainShuffle(start:min(start+miniBatchSizeErr-1, numel(idxTrainShuffle)));

            alphaBatchPIE = predictPIEOneClick(netExpert, X_valid(batchIdx,:), A_valid(batchIdx,:));
            resBatch = Y_valid(batchIdx) - alphaBatchPIE;

            dlX = dlarray(single(X_valid(batchIdx,:)'), 'CB');
            dlR = dlarray(single(resBatch'), 'CB');

            [lossErr, gradientsErr] = dlfeval(@errGradientsPIEOneClick, ...
                netErr, dlX, dlR, lambdaErrL2);

            [netErr, avgGrad, avgSqGrad] = adamupdate( ...
                netErr, gradientsErr, avgGrad, avgSqGrad, ...
                iteration, learnRateErr, gradDecay, sqGradDecay);
        end

        if mod(epoch, 250) == 0
            learnRateErr = learnRateErr * 0.5;
        end

        errTrainPred = predictErrPIEOneClick(netErr, X_valid(idxTrain,:));
        errValPred   = predictErrPIEOneClick(netErr, X_valid(idxVal,:));

        alphaTrainPIEE = clipAlphaPIEOneClick(alphaTrainPIE + errTrainPred);
        alphaValPIEE   = clipAlphaPIEOneClick(alphaValPIE   + errValPred);

        metTrainE = calcMetricsPIEOneClick_R2FIX20260513(Y_valid(idxTrain), alphaTrainPIEE);
        metValE   = calcMetricsPIEOneClick_R2FIX20260513(Y_valid(idxVal), alphaValPIEE);

        historyErr(epoch,:) = [double(gather(extractdata(lossErr))), ...
                               metValE.ARE, metValE.RMSE];

        if mod(epoch, 100) == 0 || epoch == 1 || epoch == numEpochsErr
            fprintf('[Seed %d][PIE-DF-E] Epoch %4d/%4d | Train ARE %.4f | Val ARE %.4f | Val RMSE %.4f\n', ...
                thisSeed, epoch, numEpochsErr, metTrainE.ARE, metValE.ARE, metValE.RMSE);
        end
    end

    errTrain = predictErrPIEOneClick(netErr, X_valid(idxTrain,:));
    errVal   = predictErrPIEOneClick(netErr, X_valid(idxVal,:));
    errTest  = predictErrPIEOneClick(netErr, X_valid(idxTest,:));

    alphaTrainPIEE = clipAlphaPIEOneClick(alphaTrainPIE + errTrain);
    alphaValPIEE   = clipAlphaPIEOneClick(alphaValPIE   + errVal);
    alphaTestPIEE  = clipAlphaPIEOneClick(alphaTestPIE  + errTest);

    metPIEETrain = calcMetricsPIEOneClick_R2FIX20260513(Y_valid(idxTrain), alphaTrainPIEE);
    metPIEEVal   = calcMetricsPIEOneClick_R2FIX20260513(Y_valid(idxVal), alphaValPIEE);
    metPIEETest  = calcMetricsPIEOneClick_R2FIX20260513(Y_valid(idxTest), alphaTestPIEE);

    %% --------------------------------------------------------------------
    %  记录该 seed 的测试集结果
    % ---------------------------------------------------------------------
    seedResults.PIE_DF_Test_ARE(iSeed)  = metPIETest.ARE;
    seedResults.PIE_DF_Test_MAE(iSeed)  = metPIETest.MAE;
    seedResults.PIE_DF_Test_RMSE(iSeed) = metPIETest.RMSE;
    seedResults.PIE_DF_Test_R2(iSeed)   = metPIETest.R2;

    seedResults.PIE_DF_E_Test_ARE(iSeed)  = metPIEETest.ARE;
    seedResults.PIE_DF_E_Test_MAE(iSeed)  = metPIEETest.MAE;
    seedResults.PIE_DF_E_Test_RMSE(iSeed) = metPIEETest.RMSE;
    seedResults.PIE_DF_E_Test_R2(iSeed)   = metPIEETest.R2;

    if isfinite(bestTraditionalARE) && bestTraditionalARE > 0
        seedResults.PIE_DF_ARE_Improvement_vs_BestTraditional_percent(iSeed) = ...
            100 * (bestTraditionalARE - metPIETest.ARE) / bestTraditionalARE;

        seedResults.PIE_DF_E_ARE_Improvement_vs_BestTraditional_percent(iSeed) = ...
            100 * (bestTraditionalARE - metPIEETest.ARE) / bestTraditionalARE;
    end

    fprintf('\n[Seed %d] TEST RESULTS\n', thisSeed);
    fprintf('PIE-DF   : ARE=%.6f  MAE=%.6f  RMSE=%.6f  R2=%.6f\n', ...
        metPIETest.ARE, metPIETest.MAE, metPIETest.RMSE, metPIETest.R2);
    fprintf('PIE-DF-E : ARE=%.6f  MAE=%.6f  RMSE=%.6f  R2=%.6f\n', ...
        metPIEETest.ARE, metPIEETest.MAE, metPIEETest.RMSE, metPIEETest.R2);

    % --------------------------------------------------------------------
    % 保存代表性 seed=2026 的完整状态，后面恢复给原有作图/1003点流程
    % ---------------------------------------------------------------------
    if thisSeed == representativeSeed

        haveRepresentative = true;

        netExpertRepresentative = netExpert;
        netErrRepresentative = netErr;

        historyExpertRepresentative = historyExpert;
        historyErrRepresentative = historyErr;

        alphaTrainPIERepresentative = alphaTrainPIE;
        alphaValPIERepresentative   = alphaValPIE;
        alphaTestPIERepresentative  = alphaTestPIE;

        alphaTrainPIEERepresentative = alphaTrainPIEE;
        alphaValPIEERepresentative   = alphaValPIEE;
        alphaTestPIEERepresentative  = alphaTestPIEE;

        metPIETrainRepresentative = metPIETrain;
        metPIEValRepresentative   = metPIEVal;
        metPIETestRepresentative  = metPIETest;

        metPIEETrainRepresentative = metPIEETrain;
        metPIEEValRepresentative   = metPIEEVal;
        metPIEETestRepresentative  = metPIEETest;
    end
end

% ========================================================================
% 多随机种子统计汇总：mean / std / min / max
% ========================================================================

summaryModel = ["PIE-DF"; "PIE-DF-E"];

ARE_mean = [ ...
    mean(seedResults.PIE_DF_Test_ARE, 'omitnan'); ...
    mean(seedResults.PIE_DF_E_Test_ARE, 'omitnan')];

ARE_std = [ ...
    std(seedResults.PIE_DF_Test_ARE, 0, 'omitnan'); ...
    std(seedResults.PIE_DF_E_Test_ARE, 0, 'omitnan')];

ARE_min = [ ...
    min(seedResults.PIE_DF_Test_ARE, [], 'omitnan'); ...
    min(seedResults.PIE_DF_E_Test_ARE, [], 'omitnan')];

ARE_max = [ ...
    max(seedResults.PIE_DF_Test_ARE, [], 'omitnan'); ...
    max(seedResults.PIE_DF_E_Test_ARE, [], 'omitnan')];

MAE_mean = [ ...
    mean(seedResults.PIE_DF_Test_MAE, 'omitnan'); ...
    mean(seedResults.PIE_DF_E_Test_MAE, 'omitnan')];

MAE_std = [ ...
    std(seedResults.PIE_DF_Test_MAE, 0, 'omitnan'); ...
    std(seedResults.PIE_DF_E_Test_MAE, 0, 'omitnan')];

MAE_min = [ ...
    min(seedResults.PIE_DF_Test_MAE, [], 'omitnan'); ...
    min(seedResults.PIE_DF_E_Test_MAE, [], 'omitnan')];

MAE_max = [ ...
    max(seedResults.PIE_DF_Test_MAE, [], 'omitnan'); ...
    max(seedResults.PIE_DF_E_Test_MAE, [], 'omitnan')];

RMSE_mean = [ ...
    mean(seedResults.PIE_DF_Test_RMSE, 'omitnan'); ...
    mean(seedResults.PIE_DF_E_Test_RMSE, 'omitnan')];

RMSE_std = [ ...
    std(seedResults.PIE_DF_Test_RMSE, 0, 'omitnan'); ...
    std(seedResults.PIE_DF_E_Test_RMSE, 0, 'omitnan')];

RMSE_min = [ ...
    min(seedResults.PIE_DF_Test_RMSE, [], 'omitnan'); ...
    min(seedResults.PIE_DF_E_Test_RMSE, [], 'omitnan')];

RMSE_max = [ ...
    max(seedResults.PIE_DF_Test_RMSE, [], 'omitnan'); ...
    max(seedResults.PIE_DF_E_Test_RMSE, [], 'omitnan')];

R2_mean = [ ...
    mean(seedResults.PIE_DF_Test_R2, 'omitnan'); ...
    mean(seedResults.PIE_DF_E_Test_R2, 'omitnan')];

R2_std = [ ...
    std(seedResults.PIE_DF_Test_R2, 0, 'omitnan'); ...
    std(seedResults.PIE_DF_E_Test_R2, 0, 'omitnan')];

R2_min = [ ...
    min(seedResults.PIE_DF_Test_R2, [], 'omitnan'); ...
    min(seedResults.PIE_DF_E_Test_R2, [], 'omitnan')];

R2_max = [ ...
    max(seedResults.PIE_DF_Test_R2, [], 'omitnan'); ...
    max(seedResults.PIE_DF_E_Test_R2, [], 'omitnan')];

Improve_mean = [ ...
    mean(seedResults.PIE_DF_ARE_Improvement_vs_BestTraditional_percent, 'omitnan'); ...
    mean(seedResults.PIE_DF_E_ARE_Improvement_vs_BestTraditional_percent, 'omitnan')];

Improve_std = [ ...
    std(seedResults.PIE_DF_ARE_Improvement_vs_BestTraditional_percent, 0, 'omitnan'); ...
    std(seedResults.PIE_DF_E_ARE_Improvement_vs_BestTraditional_percent, 0, 'omitnan')];

Improve_min = [ ...
    min(seedResults.PIE_DF_ARE_Improvement_vs_BestTraditional_percent, [], 'omitnan'); ...
    min(seedResults.PIE_DF_E_ARE_Improvement_vs_BestTraditional_percent, [], 'omitnan')];

Improve_max = [ ...
    max(seedResults.PIE_DF_ARE_Improvement_vs_BestTraditional_percent, [], 'omitnan'); ...
    max(seedResults.PIE_DF_E_ARE_Improvement_vs_BestTraditional_percent, [], 'omitnan')];

MultiSeed_Summary = table( ...
    summaryModel, repmat(nSeeds,2,1), ...
    ARE_mean, ARE_std, ARE_min, ARE_max, ...
    MAE_mean, MAE_std, MAE_min, MAE_max, ...
    RMSE_mean, RMSE_std, RMSE_min, RMSE_max, ...
    R2_mean, R2_std, R2_min, R2_max, ...
    Improve_mean, Improve_std, Improve_min, Improve_max, ...
    'VariableNames', { ...
    'Model','N_Seeds', ...
    'Test_ARE_Mean','Test_ARE_Std','Test_ARE_Min','Test_ARE_Max', ...
    'Test_MAE_Mean','Test_MAE_Std','Test_MAE_Min','Test_MAE_Max', ...
    'Test_RMSE_Mean','Test_RMSE_Std','Test_RMSE_Min','Test_RMSE_Max', ...
    'Test_R2_Mean','Test_R2_Std','Test_R2_Min','Test_R2_Max', ...
    'ARE_Improvement_Mean_percent','ARE_Improvement_Std_percent', ...
    'ARE_Improvement_Min_percent','ARE_Improvement_Max_percent'} );

fprintf('\n\n============================================================\n');
fprintf('多随机种子稳健性汇总（固定 train/val/test split）\n');
fprintf('============================================================\n');
disp(seedResults);
disp(MultiSeed_Summary);

fprintf('\n论文/回复审稿人最适合报告的数字：\n');
for ii = 1:height(MultiSeed_Summary)
    fprintf('%s | Test ARE = %.6f +/- %.6f [%.6f, %.6f] | ', ...
        char(MultiSeed_Summary.Model(ii)), ...
        MultiSeed_Summary.Test_ARE_Mean(ii), ...
        MultiSeed_Summary.Test_ARE_Std(ii), ...
        MultiSeed_Summary.Test_ARE_Min(ii), ...
        MultiSeed_Summary.Test_ARE_Max(ii));
    fprintf('RMSE = %.6f +/- %.6f [%.6f, %.6f] | ', ...
        MultiSeed_Summary.Test_RMSE_Mean(ii), ...
        MultiSeed_Summary.Test_RMSE_Std(ii), ...
        MultiSeed_Summary.Test_RMSE_Min(ii), ...
        MultiSeed_Summary.Test_RMSE_Max(ii));
    fprintf('R2 = %.6f +/- %.6f [%.6f, %.6f]\n', ...
        MultiSeed_Summary.Test_R2_Mean(ii), ...
        MultiSeed_Summary.Test_R2_Std(ii), ...
        MultiSeed_Summary.Test_R2_Min(ii), ...
        MultiSeed_Summary.Test_R2_Max(ii));
end

fprintf('\n参考基线：固定测试集上 ARE 最优传统模型 = %s, ARE = %.6f\n', ...
    char(bestTraditionalModelARE), bestTraditionalARE);
fprintf('若 ARE_Improvement_Min_percent > 0，则表示所有随机种子下均优于该 ARE 基线。\n');

% ========================================================================
% 保存 Excel / CSV / MAT
% ========================================================================

multiSeedFolder = 'ML_MultiSeed_Robustness';
if ~exist(multiSeedFolder, 'dir')
    mkdir(multiSeedFolder);
end

multiSeedExcel = fullfile(multiSeedFolder, 'PIE_DF_MultiSeed_Robustness.xlsx');
if exist(multiSeedExcel, 'file')
    delete(multiSeedExcel);
end

writetable(seedResults, multiSeedExcel, 'Sheet', 'PerSeed');
writetable(MultiSeed_Summary, multiSeedExcel, 'Sheet', 'Summary');

% 保存固定 split，便于复现
originalValidIndex = find(validRows);
Split = strings(Nvalid,1);
Split(idxTrain) = "Train";
Split(idxVal)   = "Validation";
Split(idxTest)  = "Test";

Fixed_Split_Table = table( ...
    (1:Nvalid)', originalValidIndex(:), groupID_valid(:), Split, ...
    'VariableNames', {'ValidSampleIndex','Original1003Index','DatabaseID','Split'} );

writetable(Fixed_Split_Table, multiSeedExcel, 'Sheet', 'FixedSplit');

MultiSeed_Config = table( ...
    2026, representativeSeed, nSeeds, string(mat2str(trainingSeeds)), ...
    bestTraditionalModelARE, bestTraditionalARE, ...
    'VariableNames', { ...
    'FixedSplitSeed','RepresentativeTrainingSeed','NTrainingSeeds','TrainingSeeds', ...
    'BestTraditionalAREModel','BestTraditionalARE'} );

writetable(MultiSeed_Config, multiSeedExcel, 'Sheet', 'Config');

writetable(seedResults, fullfile(multiSeedFolder, 'PIE_DF_MultiSeed_PerSeed.csv'));
writetable(MultiSeed_Summary, fullfile(multiSeedFolder, 'PIE_DF_MultiSeed_Summary.csv'));

save(fullfile(multiSeedFolder, 'PIE_DF_MultiSeed_Robustness.mat'), ...
    'trainingSeeds','representativeSeed','seedResults','MultiSeed_Summary', ...
    'idxTrain','idxVal','idxTest','Fixed_Split_Table','MultiSeed_Config');

fprintf('\n多随机种子结果已保存到：%s\n', multiSeedFolder);
fprintf('  - PIE_DF_MultiSeed_Robustness.xlsx\n');
fprintf('  - PIE_DF_MultiSeed_PerSeed.csv\n');
fprintf('  - PIE_DF_MultiSeed_Summary.csv\n');
fprintf('  - PIE_DF_MultiSeed_Robustness.mat\n');

% ========================================================================
% 恢复 seed=2026 的代表性网络和变量
% 这样第11节之后原有的1003点预测、权重、图和训练曲线全部仍基于2026。
% ========================================================================

if ~haveRepresentative
    error('representativeSeed=%d 不在 trainingSeeds 中，请修改设置。', representativeSeed);
end

netExpert = netExpertRepresentative;
netErr = netErrRepresentative;

historyExpert = historyExpertRepresentative;
historyErr = historyErrRepresentative;

alphaTrainPIE = alphaTrainPIERepresentative;
alphaValPIE   = alphaValPIERepresentative;
alphaTestPIE  = alphaTestPIERepresentative;

alphaTrainPIEE = alphaTrainPIEERepresentative;
alphaValPIEE   = alphaValPIEERepresentative;
alphaTestPIEE  = alphaTestPIEERepresentative;

metPIETrain = metPIETrainRepresentative;
metPIEVal   = metPIEValRepresentative;
metPIETest  = metPIETestRepresentative;

metPIEETrain = metPIEETrainRepresentative;
metPIEEVal   = metPIEEValRepresentative;
metPIEETest  = metPIEETestRepresentative;

% 保持后续代码可能引用的变量名
layersExpert = layersExpertTemplate;
layersErr = layersErrTemplate;
learnRate = learnRate0;
learnRateErr = learnRateErr0;

fprintf('\n================ 恢复代表性 seed=%d 结果，继续原有后处理 ================\n', ...
    representativeSeed);
fprintf('PIE-DF   Test: ARE=%.6f MAE=%.6f RMSE=%.6f R2=%.6f\n', ...
    metPIETest.ARE, metPIETest.MAE, metPIETest.RMSE, metPIETest.R2);
fprintf('PIE-DF-E Test: ARE=%.6f MAE=%.6f RMSE=%.6f R2=%.6f\n', ...
    metPIEETest.ARE, metPIEETest.MAE, metPIEETest.RMSE, metPIEETest.R2);

%% ========================================================================
%  11. 生成全体1003点 ML 预测
% ========================================================================

alpha_PIE_valid = predictPIEOneClick(netExpert, X_valid, A_valid);
err_valid = predictErrPIEOneClick(netErr, X_valid);
alpha_PIEE_valid = clipAlphaPIEOneClick(alpha_PIE_valid + err_valid);

alpha_PIE_all = NaN(N_total_db,1);
alpha_PIEE_all = NaN(N_total_db,1);

alpha_PIE_all(validRows) = alpha_PIE_valid;
alpha_PIEE_all(validRows) = alpha_PIEE_valid;

W_all_valid = predictWeightsPIEOneClick(netExpert, X_valid);
meanWeights = mean(W_all_valid, 1);
[meanWeightsSorted, orderW] = sort(meanWeights, 'descend');

fprintf('\n================ 全体有效样本平均权重排序 ================\n');

for ii = 1:M
    mm = orderW(ii);
    fprintf('%2d) %-22s mean weight = %.4f\n', ...
        ii, modelNames{mm}, meanWeightsSorted(ii));
end

%% ========================================================================
%  12. 只画 ML 模型在13个数据库上的图
% ========================================================================

outFolder = 'ML_Embedded_13DB_Plots';

if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

%% ========================================================================
%  11.5 PIE-DF 和 PIE-DF-E 对整体1003个数据的总误差
%       只在命令窗口输出，不保存 Excel/CSV
% ========================================================================

fprintf('\n================ PIE-DF / PIE-DF-E 整体1003点总误差（只输出窗口，不保存） ================\n');

% 强制转成列向量，避免行列方向导致计算错误
alpha_exp_all_for_relerr  = alpha_m(:);
alpha_PIE_all_for_relerr  = alpha_PIE_all(:);
alpha_PIEE_all_for_relerr = alpha_PIEE_all(:);

N_all_relerr = length(alpha_exp_all_for_relerr);
Sample_Index_relerr = (1:N_all_relerr)';

% 数据库编号
Database_ID_relerr = NaN(N_all_relerr,1);

for db_relerr = 1:length(db_counts)
    idx_db_relerr = db_start(db_relerr):db_end(db_relerr);
    Database_ID_relerr(idx_db_relerr) = db_relerr;
end

% 初始化逐点相对误差
RelErr_PIE_DF_all   = NaN(N_all_relerr,1);
RelErr_PIE_DF_E_all = NaN(N_all_relerr,1);

% 有效点筛选
valid_PIE_relerr = isfinite(alpha_exp_all_for_relerr) & isfinite(alpha_PIE_all_for_relerr) & ...
                   alpha_exp_all_for_relerr > 1e-8 & ...
                   alpha_exp_all_for_relerr >= 0 & alpha_exp_all_for_relerr <= 1 & ...
                   alpha_PIE_all_for_relerr >= 0 & alpha_PIE_all_for_relerr <= 1;

valid_PIEE_relerr = isfinite(alpha_exp_all_for_relerr) & isfinite(alpha_PIEE_all_for_relerr) & ...
                    alpha_exp_all_for_relerr > 1e-8 & ...
                    alpha_exp_all_for_relerr >= 0 & alpha_exp_all_for_relerr <= 1 & ...
                    alpha_PIEE_all_for_relerr >= 0 & alpha_PIEE_all_for_relerr <= 1;

% 逐点相对误差
RelErr_PIE_DF_all(valid_PIE_relerr) = ...
    abs((alpha_PIE_all_for_relerr(valid_PIE_relerr) - alpha_exp_all_for_relerr(valid_PIE_relerr)) ./ ...
          alpha_exp_all_for_relerr(valid_PIE_relerr));

RelErr_PIE_DF_E_all(valid_PIEE_relerr) = ...
    abs((alpha_PIEE_all_for_relerr(valid_PIEE_relerr) - alpha_exp_all_for_relerr(valid_PIEE_relerr)) ./ ...
          alpha_exp_all_for_relerr(valid_PIEE_relerr));

% 总体平均相对误差 ARE
ARE_PIE_DF_all_relerr   = mean(RelErr_PIE_DF_all(valid_PIE_relerr), 'omitnan');
ARE_PIE_DF_E_all_relerr = mean(RelErr_PIE_DF_E_all(valid_PIEE_relerr), 'omitnan');

% RMSE
RMSE_PIE_DF_all_relerr = sqrt(mean((alpha_PIE_all_for_relerr(valid_PIE_relerr) - ...
                                    alpha_exp_all_for_relerr(valid_PIE_relerr)).^2, 'omitnan'));

RMSE_PIE_DF_E_all_relerr = sqrt(mean((alpha_PIEE_all_for_relerr(valid_PIEE_relerr) - ...
                                      alpha_exp_all_for_relerr(valid_PIEE_relerr)).^2, 'omitnan'));

% MAE
MAE_PIE_DF_all_relerr = mean(abs(alpha_PIE_all_for_relerr(valid_PIE_relerr) - ...
                                  alpha_exp_all_for_relerr(valid_PIE_relerr)), 'omitnan');

MAE_PIE_DF_E_all_relerr = mean(abs(alpha_PIEE_all_for_relerr(valid_PIEE_relerr) - ...
                                    alpha_exp_all_for_relerr(valid_PIEE_relerr)), 'omitnan');

% 标准 R2
ssRes_PIE_relerr = sum((alpha_exp_all_for_relerr(valid_PIE_relerr) - ...
                        alpha_PIE_all_for_relerr(valid_PIE_relerr)).^2, 'omitnan');
ssTot_PIE_relerr = sum((alpha_exp_all_for_relerr(valid_PIE_relerr) - ...
                        mean(alpha_exp_all_for_relerr(valid_PIE_relerr), 'omitnan')).^2, 'omitnan');

if ssTot_PIE_relerr > eps
    R2_PIE_DF_all_relerr = 1 - ssRes_PIE_relerr / ssTot_PIE_relerr;
else
    R2_PIE_DF_all_relerr = NaN;
end

ssRes_PIEE_relerr = sum((alpha_exp_all_for_relerr(valid_PIEE_relerr) - ...
                         alpha_PIEE_all_for_relerr(valid_PIEE_relerr)).^2, 'omitnan');
ssTot_PIEE_relerr = sum((alpha_exp_all_for_relerr(valid_PIEE_relerr) - ...
                         mean(alpha_exp_all_for_relerr(valid_PIEE_relerr), 'omitnan')).^2, 'omitnan');

if ssTot_PIEE_relerr > eps
    R2_PIE_DF_E_all_relerr = 1 - ssRes_PIEE_relerr / ssTot_PIEE_relerr;
else
    R2_PIE_DF_E_all_relerr = NaN;
end

% Pearson corr^2，方便和 Origin 拟合 R2 对比
if sum(valid_PIE_relerr) >= 2
    R_tmp_PIE_relerr = corr(alpha_exp_all_for_relerr(valid_PIE_relerr), ...
                            alpha_PIE_all_for_relerr(valid_PIE_relerr), 'Rows', 'complete');
    R2corr_PIE_DF_all_relerr = R_tmp_PIE_relerr^2;
else
    R2corr_PIE_DF_all_relerr = NaN;
end

if sum(valid_PIEE_relerr) >= 2
    R_tmp_PIEE_relerr = corr(alpha_exp_all_for_relerr(valid_PIEE_relerr), ...
                             alpha_PIEE_all_for_relerr(valid_PIEE_relerr), 'Rows', 'complete');
    R2corr_PIE_DF_E_all_relerr = R_tmp_PIEE_relerr^2;
else
    R2corr_PIE_DF_E_all_relerr = NaN;
end

% 命令行输出
fprintf('PIE-DF   : 有效点 = %d / %d, ARE = %.5f, MAE = %.5f, RMSE = %.5f, R2 = %.5f, R2corr = %.5f\n', ...
    sum(valid_PIE_relerr), N_all_relerr, ARE_PIE_DF_all_relerr, MAE_PIE_DF_all_relerr, ...
    RMSE_PIE_DF_all_relerr, R2_PIE_DF_all_relerr, R2corr_PIE_DF_all_relerr);

fprintf('PIE-DF-E : 有效点 = %d / %d, ARE = %.5f, MAE = %.5f, RMSE = %.5f, R2 = %.5f, R2corr = %.5f\n', ...
    sum(valid_PIEE_relerr), N_all_relerr, ARE_PIE_DF_E_all_relerr, MAE_PIE_DF_E_all_relerr, ...
    RMSE_PIE_DF_E_all_relerr, R2_PIE_DF_E_all_relerr, R2corr_PIE_DF_E_all_relerr);

% 只在命令窗口输出整体1003点误差，不保存逐点误差文件
fprintf('\n================ PIE-DF / PIE-DF-E 整体1003点误差汇总 ================\n');
fprintf('总数据点数 N_total = %d\n', N_all_relerr);
fprintf('PIE-DF   整体1003点: N_valid = %d / %d, ValidRate = %.2f%%, ARE = %.6f, MAE = %.6f, RMSE = %.6f, R2 = %.6f, R2corr = %.6f\n', ...
    sum(valid_PIE_relerr), N_all_relerr, 100*sum(valid_PIE_relerr)/N_all_relerr, ...
    ARE_PIE_DF_all_relerr, MAE_PIE_DF_all_relerr, RMSE_PIE_DF_all_relerr, ...
    R2_PIE_DF_all_relerr, R2corr_PIE_DF_all_relerr);
fprintf('PIE-DF-E 整体1003点: N_valid = %d / %d, ValidRate = %.2f%%, ARE = %.6f, MAE = %.6f, RMSE = %.6f, R2 = %.6f, R2corr = %.6f\n', ...
    sum(valid_PIEE_relerr), N_all_relerr, 100*sum(valid_PIEE_relerr)/N_all_relerr, ...
    ARE_PIE_DF_E_all_relerr, MAE_PIE_DF_E_all_relerr, RMSE_PIE_DF_E_all_relerr, ...
    R2_PIE_DF_E_all_relerr, R2corr_PIE_DF_E_all_relerr);
fprintf('====================================================================\n');

axisMin = 0.0;
axisMax = 1.0;

markerSize = 22;

saveFigures = true;
showFigures = true;

fprintf('\n================ 开始绘制 ML 模型 13数据库图 ================\n');

figPIEDF_13DB_R2FINAL = plotML13DBParityPIEOneClick_R2FIX20260513( ...
    alpha_m, ...
    alpha_PIE_all, ...
    db_counts, ...
    db_start, ...
    db_end, ...
    'PIE-DF', ...
    outFolder, ...
    axisMin, ...
    axisMax, ...
    markerSize, ...
    saveFigures, ...
    showFigures);

% 再次强制从当前 PIE-DF 子图散点读取 XData/YData，覆盖标题与蓝色 R2，并重新保存图片
forceUpdate13DBFigureMetrics_R2FINAL(figPIEDF_13DB_R2FINAL, 'PIE-DF', outFolder, axisMin, axisMax);

figPIEDFE_13DB_R2FINAL = plotML13DBParityPIEOneClick_R2FIX20260513( ...
    alpha_m, ...
    alpha_PIEE_all, ...
    db_counts, ...
    db_start, ...
    db_end, ...
    'PIE-DF-E', ...
    outFolder, ...
    axisMin, ...
    axisMax, ...
    markerSize, ...
    saveFigures, ...
    showFigures);

% 同样强制修正 PIE-DF-E 子图
forceUpdate13DBFigureMetrics_R2FINAL(figPIEDFE_13DB_R2FINAL, 'PIE-DF-E', outFolder, axisMin, axisMax);


%% ========================================================================
%  12.3 保存截图这种“13个数据库预测空泡份额子图”的 Origin 数据
%      说明：本段只导出数据，不影响前面已经绘制的 13 个子图
% ========================================================================

origin13Folder = fullfile(outFolder, 'Origin_Data_13DB');
if ~exist(origin13Folder, 'dir')
    mkdir(origin13Folder);
end

axisMin_origin = axisMin;
axisMax_origin = axisMax;

% 参考线：y=x, +30%, -30%，用于 Origin 中重画每个子图的三条线
x_ref_13db = linspace(axisMin_origin, axisMax_origin, 301)';
Reference_Lines_13DB = table( ...
    x_ref_13db, ...
    x_ref_13db, ...
    1.3 .* x_ref_13db, ...
    0.7 .* x_ref_13db, ...
    'VariableNames', { ...
        'Alpha_exp_X', ...
        'Y_equal_X', ...
        'Y_plus_30_percent', ...
        'Y_minus_30_percent'} );

% ------------------------------------------------------------------------
% A. 保存 PIE-DF：对应你截图标题中的 PIE-DF 13 databases
% ------------------------------------------------------------------------
PIEDF_originFolder = fullfile(origin13Folder, 'PIE_DF');
if ~exist(PIEDF_originFolder, 'dir')
    mkdir(PIEDF_originFolder);
end

writetable(Reference_Lines_13DB, fullfile(PIEDF_originFolder, 'Reference_Lines.csv'));
writetable(Reference_Lines_13DB, fullfile(PIEDF_originFolder, 'Reference_Lines.xlsx'));

PIE_DF_excelFile = fullfile(PIEDF_originFolder, 'PIE_DF_13DB_Origin_Data.xlsx');
if exist(PIE_DF_excelFile, 'file')
    delete(PIE_DF_excelFile);
end

PIE_DF_13DB_All_Origin_Data = table();
numDB_origin = length(db_counts);

for db = 1:numDB_origin

    idx_db = (db_start(db):db_end(db))';

    alpha_exp_db = alpha_m(idx_db);
    alpha_cal_db = alpha_PIE_all(idx_db);

    % 统一转列向量，确保每个数据库 R2 与 Origin 导出数据一致
    alpha_exp_db = alpha_exp_db(:);
    alpha_cal_db = alpha_cal_db(:);

    valid_db = isfinite(alpha_exp_db) & isfinite(alpha_cal_db) & ...
               alpha_exp_db > 1e-8 & ...
               alpha_exp_db >= 0 & alpha_exp_db <= 1 & ...
               alpha_cal_db >= 0 & alpha_cal_db <= 1;

    alpha_exp_valid = alpha_exp_db(valid_db);
    alpha_cal_valid = alpha_cal_db(valid_db);
    sample_local = find(valid_db);
    sample_global = idx_db(valid_db);

    if isempty(alpha_exp_valid)
        rel_err_db = [];
        met_db.ARE = NaN;
        met_db.R2 = NaN;
        met_db.R2_corr = NaN;
        met_db.MAE = NaN;
        met_db.RMSE = NaN;
    else
        rel_err_db = abs((alpha_cal_valid - alpha_exp_valid) ./ alpha_exp_valid);
        met_db = calcMetricsPIEOneClick_R2FIX20260513(alpha_exp_valid, alpha_cal_valid);
    end

    DB_Data = table( ...
        repmat(db, numel(alpha_exp_valid), 1), ...
        sample_local(:), ...
        sample_global(:), ...
        alpha_exp_valid(:), ...
        alpha_cal_valid(:), ...
        rel_err_db(:), ...
        repmat(met_db.ARE, numel(alpha_exp_valid), 1), ...
        repmat(met_db.R2, numel(alpha_exp_valid), 1), ...
        repmat(met_db.R2_corr, numel(alpha_exp_valid), 1), ...
        repmat(met_db.MAE, numel(alpha_exp_valid), 1), ...
        repmat(met_db.RMSE, numel(alpha_exp_valid), 1), ...
        'VariableNames', { ...
            'Database_ID', ...
            'Sample_Index_In_DB', ...
            'Sample_Index_Global', ...
            'Alpha_exp', ...
            'Alpha_cal_PIE_DF', ...
            'Relative_Error', ...
            'DB_ARE', ...
            'DB_R2', ...
            'DB_R2_corr', ...
            'DB_MAE', ...
            'DB_RMSE'} );

    PIE_DF_13DB_All_Origin_Data = [PIE_DF_13DB_All_Origin_Data; DB_Data]; %#ok<AGROW>

    writetable(DB_Data, fullfile(PIEDF_originFolder, sprintf('PIE_DF_DB_%02d_Data.csv', db)));
    writetable(DB_Data, PIE_DF_excelFile, 'Sheet', sprintf('DB_%02d_Data', db));
    writetable(Reference_Lines_13DB, PIE_DF_excelFile, 'Sheet', sprintf('DB_%02d_Lines', db));
end

writetable(PIE_DF_13DB_All_Origin_Data, fullfile(PIEDF_originFolder, 'PIE_DF_13DB_All_Data.csv'));
writetable(PIE_DF_13DB_All_Origin_Data, fullfile(PIEDF_originFolder, 'PIE_DF_13DB_All_Data.xlsx'));

fprintf('\nPIE-DF 13个数据库子图 Origin 数据已保存到：\n%s\n', PIEDF_originFolder);
fprintf('总 Excel 文件：%s\n', PIE_DF_excelFile);

% ------------------------------------------------------------------------
% B. 保存 PIE-DF-E：如果不需要，可以忽略这个文件夹
% ------------------------------------------------------------------------
PIEDFE_originFolder = fullfile(origin13Folder, 'PIE_DF_E');
if ~exist(PIEDFE_originFolder, 'dir')
    mkdir(PIEDFE_originFolder);
end

writetable(Reference_Lines_13DB, fullfile(PIEDFE_originFolder, 'Reference_Lines.csv'));
writetable(Reference_Lines_13DB, fullfile(PIEDFE_originFolder, 'Reference_Lines.xlsx'));

PIE_DF_E_excelFile = fullfile(PIEDFE_originFolder, 'PIE_DF_E_13DB_Origin_Data.xlsx');
if exist(PIE_DF_E_excelFile, 'file')
    delete(PIE_DF_E_excelFile);
end

PIE_DF_E_13DB_All_Origin_Data = table();

for db = 1:numDB_origin

    idx_db = (db_start(db):db_end(db))';

    alpha_exp_db = alpha_m(idx_db);
    alpha_cal_db = alpha_PIEE_all(idx_db);

    % 统一转列向量，确保每个数据库 R2 与 Origin 导出数据一致
    alpha_exp_db = alpha_exp_db(:);
    alpha_cal_db = alpha_cal_db(:);

    valid_db = isfinite(alpha_exp_db) & isfinite(alpha_cal_db) & ...
               alpha_exp_db > 1e-8 & ...
               alpha_exp_db >= 0 & alpha_exp_db <= 1 & ...
               alpha_cal_db >= 0 & alpha_cal_db <= 1;

    alpha_exp_valid = alpha_exp_db(valid_db);
    alpha_cal_valid = alpha_cal_db(valid_db);
    sample_local = find(valid_db);
    sample_global = idx_db(valid_db);

    if isempty(alpha_exp_valid)
        rel_err_db = [];
        met_db.ARE = NaN;
        met_db.R2 = NaN;
        met_db.R2_corr = NaN;
        met_db.MAE = NaN;
        met_db.RMSE = NaN;
    else
        rel_err_db = abs((alpha_cal_valid - alpha_exp_valid) ./ alpha_exp_valid);
        met_db = calcMetricsPIEOneClick_R2FIX20260513(alpha_exp_valid, alpha_cal_valid);
    end

    DB_Data = table( ...
        repmat(db, numel(alpha_exp_valid), 1), ...
        sample_local(:), ...
        sample_global(:), ...
        alpha_exp_valid(:), ...
        alpha_cal_valid(:), ...
        rel_err_db(:), ...
        repmat(met_db.ARE, numel(alpha_exp_valid), 1), ...
        repmat(met_db.R2, numel(alpha_exp_valid), 1), ...
        repmat(met_db.R2_corr, numel(alpha_exp_valid), 1), ...
        repmat(met_db.MAE, numel(alpha_exp_valid), 1), ...
        repmat(met_db.RMSE, numel(alpha_exp_valid), 1), ...
        'VariableNames', { ...
            'Database_ID', ...
            'Sample_Index_In_DB', ...
            'Sample_Index_Global', ...
            'Alpha_exp', ...
            'Alpha_cal_PIE_DF_E', ...
            'Relative_Error', ...
            'DB_ARE', ...
            'DB_R2', ...
            'DB_R2_corr', ...
            'DB_MAE', ...
            'DB_RMSE'} );

    PIE_DF_E_13DB_All_Origin_Data = [PIE_DF_E_13DB_All_Origin_Data; DB_Data]; %#ok<AGROW>

    writetable(DB_Data, fullfile(PIEDFE_originFolder, sprintf('PIE_DF_E_DB_%02d_Data.csv', db)));
    writetable(DB_Data, PIE_DF_E_excelFile, 'Sheet', sprintf('DB_%02d_Data', db));
    writetable(Reference_Lines_13DB, PIE_DF_E_excelFile, 'Sheet', sprintf('DB_%02d_Lines', db));
end

writetable(PIE_DF_E_13DB_All_Origin_Data, fullfile(PIEDFE_originFolder, 'PIE_DF_E_13DB_All_Data.csv'));
writetable(PIE_DF_E_13DB_All_Origin_Data, fullfile(PIEDFE_originFolder, 'PIE_DF_E_13DB_All_Data.xlsx'));

fprintf('\nPIE-DF-E 13个数据库子图 Origin 数据已保存到：\n%s\n', PIEDFE_originFolder);
fprintf('总 Excel 文件：%s\n', PIE_DF_E_excelFile);

fprintf('\n================ 13数据库 Origin 数据保存完成 ================\n');

%% ========================================================================
%  13. 总体 ML 对比图
% ========================================================================

fig = figure('Color','w', 'Position',[200 120 620 540]);

if ~showFigures
    set(fig, 'Visible', 'off');
end

hold on;
box on;

validPIE = isfinite(alpha_m) & isfinite(alpha_PIE_all) & ...
           alpha_m >= 0 & alpha_m <= 1 & ...
           alpha_PIE_all >= 0 & alpha_PIE_all <= 1;

validPIEE = isfinite(alpha_m) & isfinite(alpha_PIEE_all) & ...
            alpha_m >= 0 & alpha_m <= 1 & ...
            alpha_PIEE_all >= 0 & alpha_PIEE_all <= 1;

scatter(alpha_m(validPIE), alpha_PIE_all(validPIE), 30, ...
    'MarkerEdgeColor', [0.1 0.35 0.9], ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 1.0);

scatter(alpha_m(validPIEE), alpha_PIEE_all(validPIEE), 30, ...
    'MarkerEdgeColor', [0.9 0.25 0.25], ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', 1.0);

xline = linspace(axisMin, axisMax, 200);

plot(xline, xline, 'k-', 'LineWidth', 1.3);
plot(xline, 1.3*xline, 'k--', 'LineWidth', 1.1);
plot(xline, 0.7*xline, 'k--', 'LineWidth', 1.1);

xlim([axisMin axisMax]);
ylim([axisMin axisMax]);

axis square;

xlabel('<\alpha_{exp}> , -', ...
    'FontName','Times New Roman', ...
    'FontSize', 14);

ylabel('<\alpha_{cal}> , -', ...
    'FontName','Times New Roman', ...
    'FontSize', 14);

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize', 12, ...
    'LineWidth', 1.0, ...
    'TickDir','in', ...
    'XMinorTick','on', ...
    'YMinorTick','on');

metPIE_all = calcMetricsPIEOneClick_R2FIX20260513(alpha_m(validPIE), alpha_PIE_all(validPIE));
metPIEE_all = calcMetricsPIEOneClick_R2FIX20260513(alpha_m(validPIEE), alpha_PIEE_all(validPIEE));

legend( ...
    sprintf('PIE-DF, ARE=%.3f', metPIE_all.ARE), ...
    sprintf('PIE-DF-E, ARE=%.3f', metPIEE_all.ARE), ...
    'y = x', ...
    '\pm30%', ...
    'Location','northwest');

title('ML embedded drift-flux model: all 13 databases', ...
    'FontName','Times New Roman', ...
    'FontSize', 15, ...
    'FontWeight','normal');

text(0.04, 0.96, ...
    sprintf('PIE-DF: N=%d, ARE=%.3f, RMSE=%.3f', ...
    sum(validPIE), metPIE_all.ARE, metPIE_all.RMSE), ...
    'Units','normalized', ...
    'VerticalAlignment','top', ...
    'FontName','Times New Roman', ...
    'FontSize', 11, ...
    'Color', [0.1 0.35 0.9]);

text(0.04, 0.88, ...
    sprintf('PIE-DF-E: N=%d, ARE=%.3f, RMSE=%.3f', ...
    sum(validPIEE), metPIEE_all.ARE, metPIEE_all.RMSE), ...
    'Units','normalized', ...
    'VerticalAlignment','top', ...
    'FontName','Times New Roman', ...
    'FontSize', 11, ...
    'Color', [0.9 0.25 0.25]);

exportgraphics(fig, fullfile(outFolder, 'ML_Overall_ParityPlot.png'), 'Resolution', 300);
savefig(fig, fullfile(outFolder, 'ML_Overall_ParityPlot.fig'));

%% ========================================================================
%  14. 每个数据库误差表：预分配版本，无 table warning
% ========================================================================

fprintf('\n================ 生成 ML 每个数据库误差表 ================\n');

modelList = {'PIE-DF','PIE-DF-E'};
alphaList = {alpha_PIE_all, alpha_PIEE_all};

numMLModels = length(modelList);
numDB = length(db_counts);
numRows = numMLModels * numDB;

ML_Error_Table = table( ...
    strings(numRows,1), ...
    zeros(numRows,1), ...
    zeros(numRows,1), ...
    zeros(numRows,1), ...
    zeros(numRows,1), ...
    zeros(numRows,1), ...
    zeros(numRows,1), ...
    'VariableNames', {'Model','Database','N','ARE','MAE','RMSE','R2'} );

row = 0;

for mm = 1:numMLModels

    modelNameNow = modelList{mm};
    alpha_cal_all = alphaList{mm};

    for db = 1:numDB

        idx = db_start(db):db_end(db);

        alpha_exp = alpha_m(idx);
        alpha_cal = alpha_cal_all(idx);

        valid = isfinite(alpha_exp) & isfinite(alpha_cal) & ...
                alpha_exp > 1e-8 & ...
                alpha_exp >= 0 & alpha_exp <= 1 & ...
                alpha_cal >= 0 & alpha_cal <= 1;

        met = calcMetricsPIEOneClick_R2FIX20260513(alpha_exp(valid), alpha_cal(valid));

        row = row + 1;

        ML_Error_Table.Model(row) = string(modelNameNow);
        ML_Error_Table.Database(row) = db;
        ML_Error_Table.N(row) = sum(valid);
        ML_Error_Table.ARE(row) = met.ARE;
        ML_Error_Table.MAE(row) = met.MAE;
        ML_Error_Table.RMSE(row) = met.RMSE;
        ML_Error_Table.R2(row) = met.R2;
    end
end

writetable(ML_Error_Table, fullfile(outFolder, 'ML_Error_By_Database.csv'));

fprintf('\nML 每个数据库误差表已保存：%s\n', ...
    fullfile(outFolder, 'ML_Error_By_Database.csv'));

disp(ML_Error_Table);

%% ========================================================================
%  14.5 再次在分数据库误差表后输出整体1003点总误差
%       只输出到命令窗口，不保存文件
% ========================================================================

fprintf('\n\n====================================================================\n');
fprintf('        PIE-DF / PIE-DF-E 整体1003点总误差（最终汇总）\n');
fprintf('====================================================================\n');

% 重新从最终预测向量计算一次，避免只看到分数据库表
alpha_exp_overall  = alpha_m(:);
alpha_PIE_overall  = alpha_PIE_all(:);
alpha_PIEE_overall = alpha_PIEE_all(:);
N_overall = length(alpha_exp_overall);

valid_overall_PIE = isfinite(alpha_exp_overall) & isfinite(alpha_PIE_overall) & ...
                    alpha_exp_overall > 1e-8 & ...
                    alpha_exp_overall >= 0 & alpha_exp_overall <= 1 & ...
                    alpha_PIE_overall >= 0 & alpha_PIE_overall <= 1;

valid_overall_PIEE = isfinite(alpha_exp_overall) & isfinite(alpha_PIEE_overall) & ...
                     alpha_exp_overall > 1e-8 & ...
                     alpha_exp_overall >= 0 & alpha_exp_overall <= 1 & ...
                     alpha_PIEE_overall >= 0 & alpha_PIEE_overall <= 1;

met_overall_PIE  = calcMetricsPIEOneClick_R2FIX20260513(alpha_exp_overall(valid_overall_PIE),  alpha_PIE_overall(valid_overall_PIE));
met_overall_PIEE = calcMetricsPIEOneClick_R2FIX20260513(alpha_exp_overall(valid_overall_PIEE), alpha_PIEE_overall(valid_overall_PIEE));

Overall_1003_Error_Table = table( ...
    ["PIE-DF"; "PIE-DF-E"], ...
    [N_overall; N_overall], ...
    [sum(valid_overall_PIE); sum(valid_overall_PIEE)], ...
    [100*sum(valid_overall_PIE)/N_overall; 100*sum(valid_overall_PIEE)/N_overall], ...
    [met_overall_PIE.ARE;  met_overall_PIEE.ARE], ...
    [met_overall_PIE.MAE;  met_overall_PIEE.MAE], ...
    [met_overall_PIE.RMSE; met_overall_PIEE.RMSE], ...
    [met_overall_PIE.R2;   met_overall_PIEE.R2], ...
    [met_overall_PIE.R2_corr; met_overall_PIEE.R2_corr], ...
    'VariableNames', {'Model','N_total','N_valid','ValidRate_percent','ARE_1003','MAE_1003','RMSE_1003','R2_1003','R2corr_1003'} );

disp(Overall_1003_Error_Table);

fprintf('PIE-DF   整体1003点: N_valid = %d / %d, ARE = %.6f, MAE = %.6f, RMSE = %.6f, R2 = %.6f, R2corr = %.6f\n', ...
    sum(valid_overall_PIE), N_overall, met_overall_PIE.ARE, met_overall_PIE.MAE, met_overall_PIE.RMSE, met_overall_PIE.R2, met_overall_PIE.R2_corr);
fprintf('PIE-DF-E 整体1003点: N_valid = %d / %d, ARE = %.6f, MAE = %.6f, RMSE = %.6f, R2 = %.6f, R2corr = %.6f\n', ...
    sum(valid_overall_PIEE), N_overall, met_overall_PIEE.ARE, met_overall_PIEE.MAE, met_overall_PIEE.RMSE, met_overall_PIEE.R2, met_overall_PIEE.R2_corr);
fprintf('====================================================================\n\n');

%% ========================================================================
%  15. 绘制神经网络训练过程图
% ========================================================================

trainPlotFolder = fullfile(outFolder, 'Training_Process_And_Weights');

if ~exist(trainPlotFolder, 'dir')
    mkdir(trainPlotFolder);
end

fprintf('\n================ 绘制神经网络训练过程图 ================\n');

%% 15.1 PIE-DF 权重网络训练 Loss

fig = figure('Color','w', 'Position',[200 150 650 460]);
hold on; box on;

epochExpert = (1:size(historyExpert,1))';

plot(epochExpert, historyExpert(:,1), 'LineWidth', 1.4);

xlabel('Epoch', 'FontName','Times New Roman', 'FontSize', 13);
ylabel('Loss', 'FontName','Times New Roman', 'FontSize', 13);

title('PIE-DF expert-weight network training loss', ...
    'FontName','Times New Roman', ...
    'FontSize', 15, ...
    'FontWeight','normal');

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize', 12, ...
    'LineWidth', 1.0, ...
    'TickDir','in', ...
    'XMinorTick','on', ...
    'YMinorTick','on');

exportgraphics(fig, fullfile(trainPlotFolder, 'PIE_DF_Expert_Network_Loss.png'), 'Resolution', 300);
savefig(fig, fullfile(trainPlotFolder, 'PIE_DF_Expert_Network_Loss.fig'));

%% 15.2 PIE-DF 验证集 ARE 和 RMSE

fig = figure('Color','w', 'Position',[200 150 650 460]);
hold on; box on;

yyaxis left
plot(epochExpert, historyExpert(:,3), 'LineWidth', 1.4);
ylabel('Validation ARE', 'FontName','Times New Roman', 'FontSize', 13);

yyaxis right
plot(epochExpert, historyExpert(:,4), 'LineWidth', 1.4);
ylabel('Validation RMSE', 'FontName','Times New Roman', 'FontSize', 13);

xlabel('Epoch', 'FontName','Times New Roman', 'FontSize', 13);

title('PIE-DF validation error history', ...
    'FontName','Times New Roman', ...
    'FontSize', 15, ...
    'FontWeight','normal');

legend({'ARE','RMSE'}, 'Location','northeast');

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize', 12, ...
    'LineWidth', 1.0, ...
    'TickDir','in', ...
    'XMinorTick','on', ...
    'YMinorTick','on');

exportgraphics(fig, fullfile(trainPlotFolder, 'PIE_DF_Validation_ARE_RMSE.png'), 'Resolution', 300);
savefig(fig, fullfile(trainPlotFolder, 'PIE_DF_Validation_ARE_RMSE.fig'));

%% 15.3 PIE-DF-E 残差网络训练 Loss

fig = figure('Color','w', 'Position',[200 150 650 460]);
hold on; box on;

epochErr = (1:size(historyErr,1))';

plot(epochErr, historyErr(:,1), 'LineWidth', 1.4);

xlabel('Epoch', 'FontName','Times New Roman', 'FontSize', 13);
ylabel('Loss', 'FontName','Times New Roman', 'FontSize', 13);

title('PIE-DF-E residual network training loss', ...
    'FontName','Times New Roman', ...
    'FontSize', 15, ...
    'FontWeight','normal');

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize', 12, ...
    'LineWidth', 1.0, ...
    'TickDir','in', ...
    'XMinorTick','on', ...
    'YMinorTick','on');

exportgraphics(fig, fullfile(trainPlotFolder, 'PIE_DF_E_Residual_Network_Loss.png'), 'Resolution', 300);
savefig(fig, fullfile(trainPlotFolder, 'PIE_DF_E_Residual_Network_Loss.fig'));

%% 15.4 PIE-DF-E 验证集 ARE 和 RMSE

fig = figure('Color','w', 'Position',[200 150 650 460]);
hold on; box on;

yyaxis left
plot(epochErr, historyErr(:,2), 'LineWidth', 1.4);
ylabel('Validation ARE', 'FontName','Times New Roman', 'FontSize', 13);

yyaxis right
plot(epochErr, historyErr(:,3), 'LineWidth', 1.4);
ylabel('Validation RMSE', 'FontName','Times New Roman', 'FontSize', 13);

xlabel('Epoch', 'FontName','Times New Roman', 'FontSize', 13);

title('PIE-DF-E validation error history', ...
    'FontName','Times New Roman', ...
    'FontSize', 15, ...
    'FontWeight','normal');

legend({'ARE','RMSE'}, 'Location','northeast');

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize', 12, ...
    'LineWidth', 1.0, ...
    'TickDir','in', ...
    'XMinorTick','on', ...
    'YMinorTick','on');

exportgraphics(fig, fullfile(trainPlotFolder, 'PIE_DF_E_Validation_ARE_RMSE.png'), 'Resolution', 300);
savefig(fig, fullfile(trainPlotFolder, 'PIE_DF_E_Validation_ARE_RMSE.fig'));

%% ========================================================================
%  15.5 保存图4-图7训练过程 Origin 绘图数据
%      图4: PIE-DF training loss
%      图5: PIE-DF validation ARE/RMSE
%      图6: PIE-DF-E training loss
%      图7: PIE-DF-E validation ARE/RMSE
% ========================================================================

originTrainDataFolder = fullfile(trainPlotFolder, 'Origin_Data_Fig4_to_Fig7');

if ~exist(originTrainDataFolder, 'dir')
    mkdir(originTrainDataFolder);
end

% ---------------- 图4：PIE-DF 专家权重网络训练 loss ----------------
Fig4_PIE_DF_Training_Loss_Origin = table( ...
    epochExpert(:), ...
    historyExpert(:,1), ...
    'VariableNames', { ...
        'Epoch', ...
        'PIE_DF_Training_Loss'} );

writetable(Fig4_PIE_DF_Training_Loss_Origin, ...
    fullfile(originTrainDataFolder, 'Fig4_PIE_DF_Training_Loss_Origin.csv'));

writetable(Fig4_PIE_DF_Training_Loss_Origin, ...
    fullfile(originTrainDataFolder, 'Fig4_PIE_DF_Training_Loss_Origin.xlsx'));

% ---------------- 图5：PIE-DF 验证集 ARE 和 RMSE ----------------
Fig5_PIE_DF_Validation_ARE_RMSE_Origin = table( ...
    epochExpert(:), ...
    historyExpert(:,3), ...
    historyExpert(:,4), ...
    'VariableNames', { ...
        'Epoch', ...
        'PIE_DF_Validation_ARE', ...
        'PIE_DF_Validation_RMSE'} );

writetable(Fig5_PIE_DF_Validation_ARE_RMSE_Origin, ...
    fullfile(originTrainDataFolder, 'Fig5_PIE_DF_Validation_ARE_RMSE_Origin.csv'));

writetable(Fig5_PIE_DF_Validation_ARE_RMSE_Origin, ...
    fullfile(originTrainDataFolder, 'Fig5_PIE_DF_Validation_ARE_RMSE_Origin.xlsx'));

% ---------------- 图6：PIE-DF-E 残差网络训练 loss ----------------
Fig6_PIE_DF_E_Training_Loss_Origin = table( ...
    epochErr(:), ...
    historyErr(:,1), ...
    'VariableNames', { ...
        'Epoch', ...
        'PIE_DF_E_Training_Loss'} );

writetable(Fig6_PIE_DF_E_Training_Loss_Origin, ...
    fullfile(originTrainDataFolder, 'Fig6_PIE_DF_E_Training_Loss_Origin.csv'));

writetable(Fig6_PIE_DF_E_Training_Loss_Origin, ...
    fullfile(originTrainDataFolder, 'Fig6_PIE_DF_E_Training_Loss_Origin.xlsx'));

% ---------------- 图7：PIE-DF-E 验证集 ARE 和 RMSE ----------------
Fig7_PIE_DF_E_Validation_ARE_RMSE_Origin = table( ...
    epochErr(:), ...
    historyErr(:,2), ...
    historyErr(:,3), ...
    'VariableNames', { ...
        'Epoch', ...
        'PIE_DF_E_Validation_ARE', ...
        'PIE_DF_E_Validation_RMSE'} );

writetable(Fig7_PIE_DF_E_Validation_ARE_RMSE_Origin, ...
    fullfile(originTrainDataFolder, 'Fig7_PIE_DF_E_Validation_ARE_RMSE_Origin.csv'));

writetable(Fig7_PIE_DF_E_Validation_ARE_RMSE_Origin, ...
    fullfile(originTrainDataFolder, 'Fig7_PIE_DF_E_Validation_ARE_RMSE_Origin.xlsx'));

% ---------------- 合并保存，方便 Origin 一次导入 ----------------
% historyExpert: 第1列 loss，第2列 validation MSE，第3列 validation ARE，第4列 validation RMSE
PIE_DF_Training_All_Origin = table( ...
    epochExpert(:), ...
    historyExpert(:,1), ...
    historyExpert(:,2), ...
    historyExpert(:,3), ...
    historyExpert(:,4), ...
    'VariableNames', { ...
        'Epoch', ...
        'Training_Loss', ...
        'Validation_MSE', ...
        'Validation_ARE', ...
        'Validation_RMSE'} );

% historyErr: 第1列 loss，第2列 validation ARE，第3列 validation RMSE
PIE_DF_E_Training_All_Origin = table( ...
    epochErr(:), ...
    historyErr(:,1), ...
    historyErr(:,2), ...
    historyErr(:,3), ...
    'VariableNames', { ...
        'Epoch', ...
        'Training_Loss', ...
        'Validation_ARE', ...
        'Validation_RMSE'} );

writetable(PIE_DF_Training_All_Origin, ...
    fullfile(originTrainDataFolder, 'PIE_DF_All_Training_Process_Origin.xlsx'));

writetable(PIE_DF_E_Training_All_Origin, ...
    fullfile(originTrainDataFolder, 'PIE_DF_E_All_Training_Process_Origin.xlsx'));

% 一个 Excel，多 sheet 保存图4-图7数据
combinedTrainingExcel = fullfile(originTrainDataFolder, 'Fig4_to_Fig7_All_Origin_Data.xlsx');

if exist(combinedTrainingExcel, 'file')
    delete(combinedTrainingExcel);
end

writetable(Fig4_PIE_DF_Training_Loss_Origin, ...
    combinedTrainingExcel, 'Sheet', 'Fig4_PIE_DF_Loss');

writetable(Fig5_PIE_DF_Validation_ARE_RMSE_Origin, ...
    combinedTrainingExcel, 'Sheet', 'Fig5_PIE_DF_Val');

writetable(Fig6_PIE_DF_E_Training_Loss_Origin, ...
    combinedTrainingExcel, 'Sheet', 'Fig6_PIE_DF_E_Loss');

writetable(Fig7_PIE_DF_E_Validation_ARE_RMSE_Origin, ...
    combinedTrainingExcel, 'Sheet', 'Fig7_PIE_DF_E_Val');

fprintf('\n图4-图7训练过程 Origin 数据已保存到：\n%s\n', originTrainDataFolder);
fprintf('主要文件：Fig4_to_Fig7_All_Origin_Data.xlsx\n');
fprintf('也已分别保存 Fig4/Fig5/Fig6/Fig7 的 csv 和 xlsx 文件。\n');

%% ========================================================================
%  16. 绘制总数据平均模型权重图
% ========================================================================

fprintf('\n================ 绘制总数据平均模型权重图 ================\n');

trainPlotFolder = fullfile(outFolder, 'Training_Process_And_Weights');

if ~exist(trainPlotFolder, 'dir')
    mkdir(trainPlotFolder);
end

% 计算所有有效样本的权重
W_all_valid = predictWeightsPIEOneClick(netExpert, X_valid);

% 总数据平均权重
meanWeights = mean(W_all_valid, 1, 'omitnan');

% 从大到小排序
[meanWeightsSorted, orderW] = sort(meanWeights, 'descend');

% 统一修改权重图中模型名称显示
modelNames_for_plot = modelNames;

modelNames_for_plot = strrep(modelNames_for_plot, 'Kinoshita et al.', 'Kinoshita et al.');
modelNames_for_plot = strrep(modelNames_for_plot, 'Ye et al.', 'Ye et al.');
modelNames_for_plot = strrep(modelNames_for_plot, 'Clark et al.', 'Clark et al.');
modelNames_for_plot = strrep(modelNames_for_plot, 'Chen et al.', 'Chen et al.');
modelNames_for_plot = strrep(modelNames_for_plot, 'Chexal et al.', 'Chexal et al.');
modelNames_for_plot = strrep(modelNames_for_plot, 'Maier et al.', 'Maier et al.');
modelNames_for_plot = strrep(modelNames_for_plot, 'Ozaki-Hibiki No.2', 'Ozaki-Hibiki No.2');
modelNames_for_plot = strrep(modelNames_for_plot, 'Morooka', 'Morooka');
modelNames_for_plot = strrep(modelNames_for_plot, 'Gui et al.', 'Gui et al.');
modelNames_for_plot = strrep(modelNames_for_plot, 'Paranjape et al.', 'Paranjape et al.');
modelNames_for_plot = strrep(modelNames_for_plot, 'Ren et al.', 'Ren et al.');
modelNames_for_plot = strrep(modelNames_for_plot, 'Julia et al.', 'Julia et al.');
modelNames_for_plot = strrep(modelNames_for_plot, 'Bestion', 'Bestion');

modelNames_for_plot_sorted = modelNames_for_plot(orderW);

% 作图
fig = figure('Color','w', 'Position',[100 100 950 520]);
hold on;
box on;

bar(meanWeightsSorted, 'LineWidth', 0.8);

xticks(1:M);
xticklabels(modelNames_for_plot_sorted);
xtickangle(45);

ylabel('Mean weight', ...
    'FontName','Times New Roman', ...
    'FontSize', 13);

xlabel('Traditional drift-flux model', ...
    'FontName','Times New Roman', ...
    'FontSize', 13);

title('Overall average expert weights of PIE-DF model', ...
    'FontName','Times New Roman', ...
    'FontSize', 15, ...
    'FontWeight','normal');

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize', 11, ...
    'LineWidth', 1.0, ...
    'TickDir','in', ...
    'XMinorTick','on', ...
    'YMinorTick','on');

ylim([0, max(meanWeightsSorted)*1.20 + eps]);

% 标注数值
for ii = 1:M
    text(ii, meanWeightsSorted(ii), sprintf('%.3f', meanWeightsSorted(ii)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontName','Times New Roman', ...
        'FontSize', 9);
end

% 保存图片
exportgraphics(fig, fullfile(trainPlotFolder, 'PIE_DF_Overall_Mean_Expert_Weights.png'), ...
    'Resolution', 300);

savefig(fig, fullfile(trainPlotFolder, 'PIE_DF_Overall_Mean_Expert_Weights.fig'));

% 保存权重表
Overall_Expert_Weight_Table = table( ...
    string(modelNames(orderW)), ...
    string(usedVarNames(orderW)), ...
    meanWeightsSorted(:), ...
    'VariableNames', {'Model','VariableName','MeanWeight'} );

writetable(Overall_Expert_Weight_Table, ...
    fullfile(trainPlotFolder, 'PIE_DF_Overall_Mean_Expert_Weights.csv'));

fprintf('\n总数据平均权重表：\n');
disp(Overall_Expert_Weight_Table);

fprintf('\n总数据平均权重图已保存：\n%s\n', ...
    fullfile(trainPlotFolder, 'PIE_DF_Overall_Mean_Expert_Weights.png'));
%% ========================================================================
%  17. 重新保存，包含训练曲线和权重表
% ========================================================================

save(fullfile(outFolder, 'PIE_DF_OneClick_Model_And_13DB_Predictions.mat'), ...
    'netExpert', ...
    'netErr', ...
    'modelNames', ...
    'usedVarNames', ...
    'featureNames', ...
    'muX', ...
    'stdX', ...
    'validRows', ...
    'A_valid', ...
    'X_valid', ...
    'Y_valid', ...
    'alpha_PIE_all', ...
    'alpha_PIEE_all', ...
    'alpha_m', ...
    'db_counts', ...
    'db_start', ...
    'db_end', ...
    'ML_Error_Table', ...
    'singleMetrics', ...
    'metPIETest', ...
    'metPIEETest', ...
    'historyExpert', ...
    'historyErr', ...
    'W_all_valid', ...
    'meanWeights', ...
    'Overall_Expert_Weight_Table');
fprintf('\n训练过程图和权重图已保存到文件夹：\n%s\n', trainPlotFolder);

%% ========================================================================
%  A方法最终版：只输出 PIE-DF 的 C0 / vgj 范围
%  说明：PIE-DF 按专家权重加权融合传统漂移通量模型的 C0 和 vgj。
%        PIE-DF-E 是残差网络直接修正最终空泡份额 alpha，不能唯一分解出
%        新的 C0 和 vgj，因此本版不再输出 PIE-DF-E 的 C0/vgj。
% ========================================================================

fprintf('\n================ A方法最终版：只输出 PIE-DF 的 C0 / vgj 范围 ================\n');
fprintf('说明：PIE-DF = 专家权重加权 C0/vgj；PIE-DF-E 不唯一对应 C0/vgj，本版不输出。\n');

A_paramFolder = fullfile(outFolder, 'PIE_DF_A_Method_C0_Vgj_Range');
if ~exist(A_paramFolder, 'dir')
    mkdir(A_paramFolder);
end

% 1. 准备专家权重矩阵：还原到原始 1003 点顺序
if ~exist('W_all_valid','var') || isempty(W_all_valid)
    W_all_valid = predictWeightsPIEOneClick(netExpert, X_valid);
end

if ~exist('M','var')
    M = numel(modelNames);
end

if ~exist('N_total_db','var')
    N_total_db = numel(alpha_m);
end

W_A_full = NaN(N_total_db, M);
W_A_full(validRows,:) = W_all_valid;

% 2. 按进入 PIE-DF 的 modelNames 收集对应传统模型的 C0 和 vgj
C0_A_expert  = NaN(N_total_db, M);
vgj_A_expert = NaN(N_total_db, M);
paramA_available = false(M,1);
paramA_validRate = NaN(M,1);

fprintf('\nA方法使用的专家 C0/vgj 收集情况：\n');

for mm = 1:M

    thisModel = char(modelNames{mm});
    C0_tmp = NaN(N_total_db,1);
    vgj_tmp = NaN(N_total_db,1);

    switch thisModel

        case 'Sun et al.'
            if exist('C0_Sun','var') && exist('v_gj_Sun','var')
                C0_tmp = C0_Sun(:);
                vgj_tmp = v_gj_Sun(:);
            end

        case 'Jowitt'
            if exist('C0_Jowitt','var') && exist('v_gj_Jowitt','var')
                C0_tmp = C0_Jowitt(:);
                vgj_tmp = v_gj_Jowitt(:);
            end

        case 'Bestion'
            if exist('v_gj_Bestion_all','var')
                C0_tmp = 1.2 - 0.2 .* sqrt(rho_g(:) ./ rho_f(:));
                vgj_tmp = v_gj_Bestion_all(:);
            end

        case 'Morooka'
            if exist('C0_Morooka','var') && exist('v_gj_Morooka','var')
                C0_tmp = C0_Morooka(:);
                vgj_tmp = v_gj_Morooka(:);
            end

        case 'Chexal et al.'
            if exist('C0_chexal_all','var') && exist('vgj_chexal_all','var')
                C0_tmp = C0_chexal_all(:);
                vgj_tmp = vgj_chexal_all(:);
            end

        case 'Maier et al.'
            if exist('C0_MC','var') && exist('v_gj_MC','var')
                C0_tmp = C0_MC(:);
                vgj_tmp = v_gj_MC(:);
            end

        case 'Paranjape et al.'
            C0_tmp = 1.05 .* ones(N_total_db,1);
            vgj_tmp = 0.123 .* ones(N_total_db,1);

        case 'Julia et al.'
            if exist('C0_Julia_saved','var') && exist('vgj_Julia_saved','var')
                C0_tmp = C0_Julia_saved(:);
                vgj_tmp = vgj_Julia_saved(:);
            elseif exist('C0_Julia_all','var') && exist('vgj_Julia_all','var')
                C0_tmp = C0_Julia_all(:);
                vgj_tmp = vgj_Julia_all(:);
            end

        case 'Kamei et al.'
            if exist('C0_Kamei','var') && exist('v_gj_Kamei','var')
                C0_tmp = C0_Kamei(:);
                vgj_tmp = v_gj_Kamei(:);
            end

        case 'Chen et al.'
            if exist('C0_chen_all','var') && exist('vgj_chen_all','var')
                C0_tmp = C0_chen_all(:);
                vgj_tmp = vgj_chen_all(:);
            elseif exist('C0_chen','var') && exist('v_gj_Chen','var')
                C0_tmp = C0_chen(:);
                vgj_tmp = v_gj_Chen(:);
            end

        case 'Clark et al.'
            if exist('C0_Clark_saved','var') && exist('vgj_Clark_saved','var')
                C0_tmp = C0_Clark_saved(:);
                vgj_tmp = vgj_Clark_saved(:);
            elseif exist('C0_all','var') && exist('vgj_Clark_all','var')
                C0_tmp = C0_all(:);
                vgj_tmp = vgj_Clark_all(:);
            end

        case 'Ozaki-Hibiki No.2'
            if exist('C0_0zHi2','var') && exist('v_gj_0zHi2','var')
                C0_tmp = C0_0zHi2(:);
                vgj_tmp = v_gj_0zHi2(:);
            end

        case 'Ren et al.'
            if exist('C0_ren_all','var') && exist('vgj_ren_all','var')
                C0_tmp = C0_ren_all(:);
                vgj_tmp = vgj_ren_all(:);
            end

        case 'Ye et al.'
            if exist('C0_ye_all','var') && exist('vgj_ye_all','var')
                C0_tmp = C0_ye_all(:);
                vgj_tmp = vgj_ye_all(:);
            end

        case 'Gui et al.'
            if exist('C0_gui_all','var') && exist('vgj_gui_all','var')
                C0_tmp = C0_gui_all(:);
                vgj_tmp = vgj_gui_all(:);
            end

        case 'Kinoshita et al.'
            if exist('C0_kinoshita_all','var') && exist('vgj_kinoshita_all','var')
                C0_tmp = C0_kinoshita_all(:);
                vgj_tmp = vgj_kinoshita_all(:);
            end

        case 'Hibiki-Tsukamoto'
            if exist('C0_HT','var') && exist('v_gj_HT','var')
                C0_tmp = C0_HT(:);
                vgj_tmp = v_gj_HT(:);
            end
    end

    if numel(C0_tmp) >= N_total_db && numel(vgj_tmp) >= N_total_db
        C0_tmp = C0_tmp(1:N_total_db);
        vgj_tmp = vgj_tmp(1:N_total_db);

        goodParam = isfinite(C0_tmp) & isfinite(vgj_tmp) & ...
                    C0_tmp > 0 & C0_tmp < 20 & ...
                    vgj_tmp > -50 & vgj_tmp < 100;

        paramA_validRate(mm) = mean(goodParam);

        if paramA_validRate(mm) >= 0.50
            C0_A_expert(:,mm)  = C0_tmp;
            vgj_A_expert(:,mm) = vgj_tmp;
            paramA_available(mm) = true;
            fprintf('已用于A方法：%-24s  C0/vgj有效率 = %.1f%%\n', thisModel, 100*paramA_validRate(mm));
        else
            fprintf('未用于A方法：%-24s  C0/vgj有效率过低 = %.1f%%\n', thisModel, 100*paramA_validRate(mm));
        end
    else
        fprintf('未用于A方法：%-24s  缺少完整 C0 或 vgj 向量\n', thisModel);
    end
end

% 3. 对有 C0/vgj 的专家重新归一化权重，计算每个点的等效 C0_A 和 vgj_A
C0_PIE_DF_A  = NaN(N_total_db,1);
vgj_PIE_DF_A = NaN(N_total_db,1);
N_expert_used_A = zeros(N_total_db,1);
Weight_sum_A = NaN(N_total_db,1);

for i = 1:N_total_db
    w0 = W_A_full(i,:);
    c0 = C0_A_expert(i,:);
    vg = vgj_A_expert(i,:);

    good = paramA_available(:)' & isfinite(w0) & isfinite(c0) & isfinite(vg) & w0 >= 0;

    if any(good)
        w = w0(good);
        if sum(w) > 1e-12
            w = w ./ sum(w);
            C0_PIE_DF_A(i)  = sum(w .* c0(good));
            vgj_PIE_DF_A(i) = sum(w .* vg(good));
            N_expert_used_A(i) = sum(good);
            Weight_sum_A(i) = sum(w0(good));
        end
    end
end


% 4. 用 A方法 C0/vgj 反算 alpha，检查它与 PIE-DF 主体的误差
alpha_PIE_DF_A_from_C0Vgj = NaN(N_total_db,1);
valid_alpha_A = isfinite(C0_PIE_DF_A) & isfinite(vgj_PIE_DF_A) & ...
                isfinite(j_g(:)) & isfinite(j_total(:)) & ...
                (C0_PIE_DF_A .* j_total(:) + vgj_PIE_DF_A) > 1e-12;
alpha_PIE_DF_A_from_C0Vgj(valid_alpha_A) = ...
    j_g(valid_alpha_A) ./ (C0_PIE_DF_A(valid_alpha_A) .* j_total(valid_alpha_A) + vgj_PIE_DF_A(valid_alpha_A));

met_A_alpha = calcMetricsPIEOneClick_R2FIX20260513(alpha_m(:), alpha_PIE_DF_A_from_C0Vgj(:));
fprintf('\nA方法 C0/vgj 反算 PIE-DF alpha 的整体检查：ARE=%.5f, MAE=%.5f, RMSE=%.5f, R2=%.5f\n', ...
    met_A_alpha.ARE, met_A_alpha.MAE, met_A_alpha.RMSE, met_A_alpha.R2);

% 5. 最终输出：PIE-DF_A 的 C0/vgj 范围，并把两相表观速度、空泡份额、压力放在一起
%    说明：1行4列图只画 j_f、j_g、alpha_exp、P；共用同一个 y 轴标签，只有第一个子图显示作者/数据库标签。
numDB_A = length(db_counts);
PIE_DF_A_C0_Vgj_Range_By_DB_FINAL = table();

% 压力统一换算为 MPa。你的前面代码里 P = data(:,5)*1e6，是 Pa；
% 如果后续版本中 P 已经是 MPa，这里会自动识别。
if exist('P','var')
    P_tmp_A = P(:);
    P_abs_med_A = median(abs(P_tmp_A(isfinite(P_tmp_A))), 'omitnan');
    if isfinite(P_abs_med_A) && P_abs_med_A > 1e5
        P_MPa_A = P_tmp_A ./ 1e6;
    else
        P_MPa_A = P_tmp_A;
    end
elseif exist('data','var') && size(data,2) >= 5
    P_MPa_A = data(:,5);
else
    P_MPa_A = NaN(N_total_db,1);
end

P_MPa_A = P_MPa_A(:);
if numel(P_MPa_A) ~= N_total_db
    warning('P_MPa_A 长度与总数据点数不一致，压力统计将用 NaN 填充。');
    P_MPa_A = NaN(N_total_db,1);
end

% 两个速度：液相表观速度 j_f 和气相表观速度 j_g
Jf_A = j_f(:);
Jg_A = j_g(:);
Alpha_exp_A = alpha_m(:);

if numel(Jf_A) ~= N_total_db
    warning('j_f 长度与总数据点数不一致，j_f 统计将用 NaN 填充。');
    Jf_A = NaN(N_total_db,1);
end
if numel(Jg_A) ~= N_total_db
    warning('j_g 长度与总数据点数不一致，j_g 统计将用 NaN 填充。');
    Jg_A = NaN(N_total_db,1);
end

fprintf('\n================ 最终输出：PIE-DF 每个数据库 C0 / vgj 范围，并统计 j_f / j_g / alpha / P ================\n');
fprintf('说明：1行4列图为 j_f、j_g、alpha_exp、Pressure，共用 y 轴，作者/数据库标签只放第一个子图。\n');

for db = 1:numDB_A
    idx_db = db_start(db):db_end(db);

    c0_pie = C0_PIE_DF_A(idx_db);
    vg_pie = vgj_PIE_DF_A(idx_db);
    jf_db = Jf_A(idx_db);
    jg_db = Jg_A(idx_db);
    alpha_exp_db_A = Alpha_exp_A(idx_db);
    p_db = P_MPa_A(idx_db);

    valid_pie = isfinite(c0_pie) & isfinite(vg_pie);
    valid_jf = isfinite(jf_db);
    valid_jg = isfinite(jg_db);
    valid_alpha_exp_A = isfinite(alpha_exp_db_A);
    valid_p = isfinite(p_db);

    if any(valid_pie)
        C0_min_pie = min(c0_pie(valid_pie));
        C0_max_pie = max(c0_pie(valid_pie));
        C0_mean_pie = mean(c0_pie(valid_pie), 'omitnan');
        C0_std_pie = std(c0_pie(valid_pie), 'omitnan');
        Vgj_min_pie = min(vg_pie(valid_pie));
        Vgj_max_pie = max(vg_pie(valid_pie));
        Vgj_mean_pie = mean(vg_pie(valid_pie), 'omitnan');
        Vgj_std_pie = std(vg_pie(valid_pie), 'omitnan');
    else
        C0_min_pie = NaN; C0_max_pie = NaN; C0_mean_pie = NaN; C0_std_pie = NaN;
        Vgj_min_pie = NaN; Vgj_max_pie = NaN; Vgj_mean_pie = NaN; Vgj_std_pie = NaN;
    end

    if any(valid_jf)
        Jf_min = min(jf_db(valid_jf));
        Jf_max = max(jf_db(valid_jf));
        Jf_mean = mean(jf_db(valid_jf), 'omitnan');
        Jf_std = std(jf_db(valid_jf), 'omitnan');
    else
        Jf_min = NaN; Jf_max = NaN; Jf_mean = NaN; Jf_std = NaN;
    end

    if any(valid_jg)
        Jg_min = min(jg_db(valid_jg));
        Jg_max = max(jg_db(valid_jg));
        Jg_mean = mean(jg_db(valid_jg), 'omitnan');
        Jg_std = std(jg_db(valid_jg), 'omitnan');
    else
        Jg_min = NaN; Jg_max = NaN; Jg_mean = NaN; Jg_std = NaN;
    end

    if any(valid_alpha_exp_A)
        Alpha_exp_min = min(alpha_exp_db_A(valid_alpha_exp_A));
        Alpha_exp_max = max(alpha_exp_db_A(valid_alpha_exp_A));
        Alpha_exp_mean = mean(alpha_exp_db_A(valid_alpha_exp_A), 'omitnan');
        Alpha_exp_std = std(alpha_exp_db_A(valid_alpha_exp_A), 'omitnan');
    else
        Alpha_exp_min = NaN; Alpha_exp_max = NaN; Alpha_exp_mean = NaN; Alpha_exp_std = NaN;
    end

    if any(valid_p)
        P_min_MPa = min(p_db(valid_p));
        P_max_MPa = max(p_db(valid_p));
        P_mean_MPa = mean(p_db(valid_p), 'omitnan');
        P_std_MPa = std(p_db(valid_p), 'omitnan');
    else
        P_min_MPa = NaN; P_max_MPa = NaN; P_mean_MPa = NaN; P_std_MPa = NaN;
    end

    fprintf('DB %02d | PIE-DF_A | N=%d, N_valid=%d | C0=[%.4f, %.4f], mean=%.4f | vgj=[%.4f, %.4f], mean=%.4f | jf=[%.4f, %.4f] | jg=[%.4f, %.4f] | alpha=[%.4f, %.4f] | P=[%.4f, %.4f] MPa\n', ...
        db, db_counts(db), sum(valid_pie), ...
        C0_min_pie, C0_max_pie, C0_mean_pie, ...
        Vgj_min_pie, Vgj_max_pie, Vgj_mean_pie, ...
        Jf_min, Jf_max, Jg_min, Jg_max, Alpha_exp_min, Alpha_exp_max, P_min_MPa, P_max_MPa);

    PIE_DF_A_C0_Vgj_Range_By_DB_FINAL = [PIE_DF_A_C0_Vgj_Range_By_DB_FINAL; table( ...
        string('PIE-DF_A'), db, db_counts(db), sum(valid_pie), ...
        Jf_min, Jf_max, Jf_mean, Jf_std, ...
        Jg_min, Jg_max, Jg_mean, Jg_std, ...
        Alpha_exp_min, Alpha_exp_max, Alpha_exp_mean, Alpha_exp_std, ...
        P_min_MPa, P_max_MPa, P_mean_MPa, P_std_MPa, ...
        C0_min_pie, C0_max_pie, C0_mean_pie, C0_std_pie, ...
        Vgj_min_pie, Vgj_max_pie, Vgj_mean_pie, Vgj_std_pie, ...
        'VariableNames', { ...
            'Model', 'Database_ID', 'N_total', 'N_valid_C0_Vgj', ...
            'Jf_min', 'Jf_max', 'Jf_mean', 'Jf_std', ...
            'Jg_min', 'Jg_max', 'Jg_mean', 'Jg_std', ...
            'Alpha_exp_min', 'Alpha_exp_max', 'Alpha_exp_mean', 'Alpha_exp_std', ...
            'Pressure_MPa_min', 'Pressure_MPa_max', 'Pressure_MPa_mean', 'Pressure_MPa_std', ...
            'C0_min', 'C0_max', 'C0_mean', 'C0_std', ...
            'Vgj_min', 'Vgj_max', 'Vgj_mean', 'Vgj_std'} )];
end

% 6. 每个点的 PIE-DF C0/vgj 明细，同时包含 j_f、j_g、alpha、P
PointIndex_A = (1:N_total_db)';
Database_ID_A = NaN(N_total_db,1);
for db = 1:numDB_A
    Database_ID_A(db_start(db):db_end(db)) = db;
end

PIE_DF_A_C0_Vgj_Pointwise_FINAL = table( ...
    PointIndex_A, ...
    Database_ID_A, ...
    Jf_A(:), ...
    Jg_A(:), ...
    Alpha_exp_A(:), ...
    P_MPa_A(:), ...
    alpha_PIE_all(:), ...
    C0_PIE_DF_A(:), ...
    vgj_PIE_DF_A(:), ...
    alpha_PIE_DF_A_from_C0Vgj(:), ...
    N_expert_used_A(:), ...
    Weight_sum_A(:), ...
    'VariableNames', { ...
        'PointIndex', ...
        'Database_ID', ...
        'Jf_m_per_s', ...
        'Jg_m_per_s', ...
        'Alpha_exp', ...
        'Pressure_MPa', ...
        'Alpha_PIE_DF', ...
        'C0_PIE_DF_A', ...
        'Vgj_PIE_DF_A', ...
        'Alpha_from_A_C0_Vgj', ...
        'N_expert_used_A', ...
        'Weight_sum_used_A'} );

% 7. 输出专家参数可用性表
Expert_C0_Vgj_Availability_FINAL = table( ...
    string(modelNames(:)), ...
    paramA_available(:), ...
    paramA_validRate(:), ...
    'VariableNames', {'Model', 'Used_in_A_Method', 'C0_Vgj_Valid_Rate'} );

% 8. 生成 1 行 4 列范围图：j_g、j_f、alpha_exp、Pressure
%    风格与原来的三联图一致：彩色水平范围条 + 黑色均值点 + 共用 y 轴作者标签。
try
    % 作者 / 数据库标签：只放在第一个子图，后面三个子图共用同一 y 轴但不重复显示
    if exist('dbAuthorNames','var') && numel(dbAuthorNames) >= numDB_A
        baseAuthorLabels_A = cellstr(string(dbAuthorNames(:)));
        baseAuthorLabels_A = baseAuthorLabels_A(1:numDB_A);
    else
        % 如果原程序里没有 dbAuthorNames，就使用与论文图一致的默认标签
        defaultAuthors_A = { ...
            'Anklam et al.'; ...
            'Kumamaru et al.'; ...
            'Morooka et al.'; ...
            'Morooka et al.'; ...
            'Inoue et al.'; ...
            'Gui et al.'; ...
            'Yun et al.'; ...
            'Griffiths'; ...
            'Chen et al.'; ...
            'Kamei et al.'; ...
            'Kamei et al.'; ...
            'Shen et al.'; ...
            'Yang et al.'};
        if numel(defaultAuthors_A) >= numDB_A
            baseAuthorLabels_A = defaultAuthors_A(1:numDB_A);
        else
            baseAuthorLabels_A = arrayfun(@(x) sprintf('DB %02d', x), 1:numDB_A, 'UniformOutput', false);
        end
    end

    yLabels_A = cell(numDB_A,1);
    for dbp = 1:numDB_A
        yLabels_A{dbp} = sprintf('%s,  N = %d', baseAuthorLabels_A{dbp}, db_counts(dbp));
    end

    % y 轴位置递增，配合 YDir reverse，使第 1 个数据库在最上面
    yPos_A = 1:numDB_A;

    % 13 个数据库使用 13 种固定颜色，避免 MATLAB 默认 7 色循环造成重复
    barColors_A_full = [ ...
        0.1216 0.4667 0.7059;  ... % 01 blue
        1.0000 0.4980 0.0549;  ... % 02 orange
        0.1725 0.6275 0.1725;  ... % 03 green
        0.8392 0.1529 0.1569;  ... % 04 red
        0.5804 0.4039 0.7412;  ... % 05 purple
        0.5490 0.3373 0.2941;  ... % 06 brown
        0.8902 0.4667 0.7608;  ... % 07 pink
        0.4980 0.4980 0.4980;  ... % 08 gray
        0.7373 0.7412 0.1333;  ... % 09 olive
        0.0902 0.7451 0.8118;  ... % 10 cyan
        0.0000 0.2471 0.4549;  ... % 11 dark blue
        0.7686 0.3059 0.0000;  ... % 12 dark orange
        0.0000 0.3922 0.0000];     % 13 dark green
    if numDB_A <= size(barColors_A_full,1)
        barColors_A = barColors_A_full(1:numDB_A,:);
    else
        barColors_A = lines(numDB_A);
    end

    figRange_A = figure('Color','w', 'Position', [60 80 1600 760]);
    tl_A = tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact');

    % 顺序按你参考图风格：气相速度、液相速度、空泡份额、压力
    panelTitles_A = { ...
        '(a) Gas superficial velocity', ...
        '(b) Liquid superficial velocity', ...
        '(c) Void fraction', ...
        '(d) Pressure'};

    % j_f 和 j_g 使用 TeX 斜体；alpha 也使用 TeX 形式
    xLabels_A = { ...
        '{\it j}_g (m/s)', ...
        '{\it j}_f (m/s)', ...
        '<\alpha> (-)', ...
        'Pressure (MPa)'};

    minCols_A = {'Jg_min', 'Jf_min', 'Alpha_exp_min', 'Pressure_MPa_min'};
    maxCols_A = {'Jg_max', 'Jf_max', 'Alpha_exp_max', 'Pressure_MPa_max'};
    meanCols_A = {'Jg_mean', 'Jf_mean', 'Alpha_exp_mean', 'Pressure_MPa_mean'};

    for ip = 1:4
        axp = nexttile(tl_A, ip);
        hold(axp, 'on');
        box(axp, 'on');

        xMin_A = PIE_DF_A_C0_Vgj_Range_By_DB_FINAL.(minCols_A{ip});
        xMax_A = PIE_DF_A_C0_Vgj_Range_By_DB_FINAL.(maxCols_A{ip});
        xMean_A = PIE_DF_A_C0_Vgj_Range_By_DB_FINAL.(meanCols_A{ip});

        for dbp = 1:numDB_A
            yy = yPos_A(dbp);
            x1 = xMin_A(dbp);
            x2 = xMax_A(dbp);
            xm = xMean_A(dbp);

            if isfinite(x1) && isfinite(x2)
                if x2 < x1
                    tmpx = x1; x1 = x2; x2 = tmpx;
                end

                % 彩色水平范围条，和你原来的风格一致
                hBar = 0.36;
                patch(axp, ...
                    [x1 x2 x2 x1], ...
                    [yy-hBar yy-hBar yy+hBar yy+hBar], ...
                    barColors_A(dbp,:), ...
                    'EdgeColor','k', ...
                    'LineWidth',0.7, ...
                    'FaceAlpha',0.95);
            end

            % 黑色均值点
            if isfinite(xm)
                plot(axp, xm, yy, 'ko', ...
                    'MarkerFaceColor','k', ...
                    'MarkerSize',4.2, ...
                    'LineWidth',0.7);
            end
        end

        ylim(axp, [0.4, numDB_A+0.6]);
        yticks(axp, 1:numDB_A);
        set(axp, 'YDir','reverse');

        if ip == 1
            yticklabels(axp, yLabels_A);
        else
            yticklabels(axp, repmat({''}, numDB_A, 1));
        end

        xlabel(axp, xLabels_A{ip}, 'FontName','Times New Roman', 'FontSize', 14, 'Interpreter','tex');
        title(axp, panelTitles_A{ip}, ...
            'FontName','Times New Roman', ...
            'FontSize', 15, ...
            'FontWeight','normal');

        set(axp, ...
            'FontName','Times New Roman', ...
            'FontSize', 12, ...
            'LineWidth', 1.0, ...
            'TickDir','in', ...
            'XGrid','on', ...
            'YGrid','on', ...
            'XMinorGrid','on', ...
            'YMinorGrid','on', ...
            'GridLineStyle','-', ...
            'MinorGridLineStyle',':', ...
            'Layer','top');

        % 给 x 轴留一点边距，压力/速度等范围看起来更像原图
        finiteX_A = [xMin_A(:); xMax_A(:); xMean_A(:)];
        finiteX_A = finiteX_A(isfinite(finiteX_A));
        if ~isempty(finiteX_A)
            xx1 = min(finiteX_A);
            xx2 = max(finiteX_A);
            if abs(xx2 - xx1) < eps
                padX = max(abs(xx2),1) * 0.05;
            else
                padX = 0.04 * (xx2 - xx1);
            end
            xlim(axp, [max(0, xx1-padX), xx2+padX]);
        end
    end

    sgtitle(figRange_A, 'Ranges of {\it j}_g, {\it j}_f, void fraction and pressure in 13 databases', ...
        'FontName','Times New Roman', 'FontSize', 17, 'FontWeight','normal', 'Interpreter','tex');

    exportgraphics(figRange_A, fullfile(A_paramFolder, 'FINAL_DB_Ranges_jg_jf_alpha_pressure_1x4_BARSTYLE_ITALIC_13COLORS.png'), 'Resolution', 300);
    savefig(figRange_A, fullfile(A_paramFolder, 'FINAL_DB_Ranges_jg_jf_alpha_pressure_1x4_BARSTYLE_ITALIC_13COLORS.fig'));
catch ME_plot_A
    warning('1行4列 jg/jf/alpha/pressure 彩色范围图生成失败：%s', ME_plot_A.message);
end

% 9. 保存最终版结果：PIE-DF C0/vgj + 1行4列图用的 jf/jg/alpha/P 范围
writetable(PIE_DF_A_C0_Vgj_Range_By_DB_FINAL, ...
    fullfile(A_paramFolder, 'FINAL_PIE_DF_ONLY_C0_Vgj_and_jf_jg_alpha_P_Range_By_Database.csv'));
writetable(PIE_DF_A_C0_Vgj_Range_By_DB_FINAL, ...
    fullfile(A_paramFolder, 'FINAL_PIE_DF_ONLY_C0_Vgj_and_jf_jg_alpha_P_Range_By_Database.xlsx'));

writetable(PIE_DF_A_C0_Vgj_Pointwise_FINAL, ...
    fullfile(A_paramFolder, 'FINAL_PIE_DF_ONLY_C0_Vgj_and_jf_jg_alpha_P_Pointwise.csv'));
writetable(PIE_DF_A_C0_Vgj_Pointwise_FINAL, ...
    fullfile(A_paramFolder, 'FINAL_PIE_DF_ONLY_C0_Vgj_and_jf_jg_alpha_P_Pointwise.xlsx'));

writetable(Expert_C0_Vgj_Availability_FINAL, ...
    fullfile(A_paramFolder, 'FINAL_PIE_DF_ONLY_Expert_C0_Vgj_Availability.csv'));
writetable(Expert_C0_Vgj_Availability_FINAL, ...
    fullfile(A_paramFolder, 'FINAL_PIE_DF_ONLY_Expert_C0_Vgj_Availability.xlsx'));

save(fullfile(A_paramFolder, 'FINAL_PIE_DF_ONLY_A_Method_C0_Vgj_Result.mat'), ...
    'C0_A_expert', 'vgj_A_expert', 'paramA_available', 'paramA_validRate', ...
    'W_A_full', 'C0_PIE_DF_A', 'vgj_PIE_DF_A', ...
    'PIE_DF_A_C0_Vgj_Range_By_DB_FINAL', 'PIE_DF_A_C0_Vgj_Pointwise_FINAL', ...
    'Expert_C0_Vgj_Availability_FINAL');

fprintf('\n最终版 PIE-DF-only 结果已保存到：\n%s\n', A_paramFolder);
fprintf('主要文件：\n');
fprintf('  1) FINAL_PIE_DF_ONLY_C0_Vgj_and_jf_jg_alpha_P_Range_By_Database.xlsx\n');
fprintf('  2) FINAL_PIE_DF_ONLY_C0_Vgj_and_jf_jg_alpha_P_Pointwise.xlsx\n');
fprintf('  3) FINAL_PIE_DF_ONLY_Expert_C0_Vgj_Availability.xlsx\n');
fprintf('  4) FINAL_DB_Ranges_jf_jg_alpha_pressure_1x4.png / .fig\n');
fprintf('====================================================================\n');

%% ========================================================================
%                           Local functions
% ========================================================================

function [loss, gradients] = expertGradientsPIEOneClick(net, dlX, dlA, dlY, dlWtrue, lambdaCE, lambdaL2W)

    logits = forward(net, dlX);

    W = softmaxColumnsPIEOneClick(logits);

    alphaHat = sum(W .* dlA, 1);

    lossPred = mean((alphaHat - dlY).^2, 'all');

    lossCE = -mean(sum(dlWtrue .* log(W + 1e-8), 1), 'all');

    lossL2W = mean(W.^2, 'all');

    loss = lossPred + lambdaCE * lossCE + lambdaL2W * lossL2W;

    gradients = dlgradient(loss, net.Learnables);
end

function [loss, gradients] = errGradientsPIEOneClick(net, dlX, dlR, lambdaErrL2)

    errHat = forward(net, dlX);

    lossMSE = mean((errHat - dlR).^2, 'all');

    lossL2 = mean(errHat.^2, 'all');

    loss = lossMSE + lambdaErrL2 * lossL2;

    gradients = dlgradient(loss, net.Learnables);
end

function W = predictWeightsPIEOneClick(net, X)

    dlX = dlarray(single(X'), 'CB');

    logits = predict(net, dlX);

    Wdl = softmaxColumnsPIEOneClick(logits);

    W = gather(extractdata(Wdl))';
end

function alphaHat = predictPIEOneClick(net, X, A)

    W = predictWeightsPIEOneClick(net, X);

    alphaHat = sum(W .* A, 2);

    alphaHat = clipAlphaPIEOneClick(alphaHat);
end

function err = predictErrPIEOneClick(net, X)

    dlX = dlarray(single(X'), 'CB');

    dlErr = predict(net, dlX);

    err = gather(extractdata(dlErr))';

    err = err(:);

    err = max(min(err, 0.25), -0.25);
end

function W = softmaxColumnsPIEOneClick(Z)

    Z = Z - max(Z, [], 1);

    EZ = exp(Z);

    W = EZ ./ sum(EZ, 1);
end

function y = clipAlphaPIEOneClick(y)

    y = min(max(y(:), 1e-5), 0.999);
end

function met = calcMetricsPIEOneClick_R2FIX20260513(yTrue, yPred)

    % 统一转为列向量，避免行/列向量隐式扩展导致 R2、ARE、RMSE 算错
    yTrue = yTrue(:);
    yPred = yPred(:);

    % 基础有效性筛选：真实值必须大于 0，避免 ARE 分母为 0
    valid = isfinite(yTrue) & isfinite(yPred) & yTrue > 1e-8;

    yTrue = yTrue(valid);
    yPred = yPred(valid);

    if isempty(yTrue)
        met.MSE = NaN;
        met.ARE = NaN;
        met.MAE = NaN;
        met.RMSE = NaN;
        met.R2 = NaN;
        met.R2_corr = NaN;
        return;
    end

    err = yPred - yTrue;

    met.MSE  = mean(err.^2, 'omitnan');
    met.ARE  = mean(abs(err ./ yTrue), 'omitnan');
    met.MAE  = mean(abs(err), 'omitnan');
    met.RMSE = sqrt(met.MSE);

    % 标准决定系数：R2 = 1 - SS_res / SS_tot
    % 注意：这里 yTrue/yPred 已经强制为列向量，避免 MATLAB 隐式扩展成矩阵
    ssRes = sum((yTrue - yPred).^2, 'omitnan');
    ssTot = sum((yTrue - mean(yTrue, 'omitnan')).^2, 'omitnan');

    if ssTot <= eps
        met.R2 = NaN;
    else
        met.R2 = 1 - ssRes / ssTot;
    end

    % 相关系数平方：仅用于和 Origin 线性拟合 R2 或论文中 corr^2 对比
    if numel(yTrue) >= 2
        R_tmp = corr(yTrue, yPred, 'Rows', 'complete');
        met.R2_corr = R_tmp.^2;
    else
        met.R2_corr = NaN;
    end
end

function fig = plotML13DBParityPIEOneClick_R2FIX20260513(alpha_exp_all, alpha_cal_all, db_counts, db_start, db_end, ...
    modelName, outFolder, axisMin, axisMax, markerSize, saveFigures, showFigures)

    % 终极修正版：R2/ARE/RMSE 只从当前子图真正画出来的散点 XData/YData 计算。
    % 这样无论 alpha_exp_all / alpha_cal_all 原来是行向量还是列向量，
    % 也无论前面有没有旧变量，图上数字都和可见散点完全一致。

    fig = figure('Color','w', 'Position', [80 50 1450 900]);

    if ~showFigures
        set(fig, 'Visible', 'off');
    end

    tiledlayout(3,5, 'TileSpacing','compact', 'Padding','compact');

    for db = 1:length(db_counts)

        nexttile;

        idx = db_start(db):db_end(db);

        alpha_exp_raw = alpha_exp_all(idx);
        alpha_cal_raw = alpha_cal_all(idx);

        alpha_exp_raw = alpha_exp_raw(:);
        alpha_cal_raw = alpha_cal_raw(:);

        valid = isfinite(alpha_exp_raw) & isfinite(alpha_cal_raw) & ...
                alpha_exp_raw >= 0 & alpha_exp_raw <= 1 & ...
                alpha_cal_raw >= 0 & alpha_cal_raw <= 1;

        alpha_exp_plot = alpha_exp_raw(valid);
        alpha_cal_plot = alpha_cal_raw(valid);

        hold on;
        box on;

        % 先画散点，并从散点对象中重新读取 XData/YData 来计算指标
        hScat = scatter(alpha_exp_plot, alpha_cal_plot, markerSize, ...
            'MarkerEdgeColor', [1.0 0.55 0.75], ...
            'MarkerFaceColor', 'none', ...
            'LineWidth', 0.9);

        x_plot = hScat.XData(:);
        y_plot = hScat.YData(:);

        validPlot = isfinite(x_plot) & isfinite(y_plot) & x_plot > 1e-8;
        x_plot = x_plot(validPlot);
        y_plot = y_plot(validPlot);

        if isempty(x_plot)
            ARE_plot = NaN;
            RMSE_plot = NaN;
            R2_plot = NaN;
            R2_corr_plot = NaN;
        else
            ARE_plot = mean(abs((y_plot - x_plot) ./ x_plot), 'omitnan');
            RMSE_plot = sqrt(mean((y_plot - x_plot).^2, 'omitnan'));

            ssRes = sum((x_plot - y_plot).^2, 'omitnan');
            ssTot = sum((x_plot - mean(x_plot, 'omitnan')).^2, 'omitnan');

            if ssTot > eps
                R2_plot = 1 - ssRes / ssTot;
            else
                R2_plot = NaN;
            end

            if numel(x_plot) >= 2
                R_tmp = corr(x_plot, y_plot, 'Rows', 'complete');
                R2_corr_plot = R_tmp.^2;
            else
                R2_corr_plot = NaN;
            end
        end

        fprintf('FORCE_PLOT_R2 | %s DB %02d: N=%d, ARE=%.5f, RMSE=%.5f, R2=%.5f, R2_corr=%.5f\n', ...
            modelName, db, numel(x_plot), ARE_plot, RMSE_plot, R2_plot, R2_corr_plot);

        xline = linspace(axisMin, axisMax, 200);

        plot(xline, xline, 'b-', 'LineWidth', 1.1);
        plot(xline, 1.3*xline, 'r--', 'LineWidth', 1.1);
        plot(xline, 0.7*xline, 'r--', 'LineWidth', 1.1);

        xlim([axisMin axisMax]);
        ylim([axisMin axisMax]);

        axis square;

        title(sprintf('DB %d, N=%d\nARE=%.3f, R^2=%.3f', ...
            db, numel(x_plot), ARE_plot, R2_plot), ...
            'FontName','Times New Roman', ...
            'FontSize', 9.5, ...
            'FontWeight','normal');

        set(gca, ...
            'FontName','Times New Roman', ...
            'FontSize', 8.5, ...
            'LineWidth', 0.8, ...
            'TickDir','in', ...
            'XMinorTick','on', ...
            'YMinorTick','on');

        if db > 10
            xlabel('\alpha_{exp}', 'FontName','Times New Roman');
        end

        if mod(db-1,5) == 0
            ylabel('\alpha_{cal}', 'FontName','Times New Roman');
        end

        text(0.08, 0.88, '+30%', ...
            'Units','normalized', ...
            'FontName','Times New Roman', ...
            'FontSize', 8, ...
            'Color','k');

        text(0.66, 0.30, '-30%', ...
            'Units','normalized', ...
            'FontName','Times New Roman', ...
            'FontSize', 8, ...
            'Color','k');

        text(0.05, 0.05, sprintf('R^2 = %.3f', R2_plot), ...
            'Units','normalized', ...
            'FontName','Times New Roman', ...
            'FontSize', 8.5, ...
            'Color','b');
    end

    sgtitle(sprintf('%s: ML embedded drift-flux model in 13 databases  [R2 FINAL]', modelName), ...
        'FontName','Times New Roman', ...
        'FontSize', 17, ...
        'FontWeight','normal');

    if saveFigures

        safeName = makeSafeFileNamePIEOneClick(modelName);

        exportgraphics(fig, fullfile(outFolder, sprintf('%s_13DB_Overview.png', safeName)), ...
            'Resolution', 300);

        savefig(fig, fullfile(outFolder, sprintf('%s_13DB_Overview.fig', safeName)));
    end
end

% 保留原函数名的包装器：如果主程序或你的旧代码还在调用旧名字，也会使用同一套强制计算逻辑。
function plotML13DBParityPIEOneClick(alpha_exp_all, alpha_cal_all, db_counts, db_start, db_end, ...
    modelName, outFolder, axisMin, axisMax, markerSize, saveFigures, showFigures)

    plotML13DBParityPIEOneClick_R2FIX20260513(alpha_exp_all, alpha_cal_all, db_counts, db_start, db_end, ...
        modelName, outFolder, axisMin, axisMax, markerSize, saveFigures, showFigures);
end

function forceUpdate13DBFigureMetrics_R2FINAL(fig, modelName, outFolder, axisMin, axisMax)

    % 这个函数只读取 figure 里面已经画出来的 scatter 点，强制覆盖每个子图的 ARE 和 R2。
    % 图上看到什么点，就按什么点算，避免任何旧变量/旧函数/行列向量问题。

    if isempty(fig) || ~ishandle(fig)
        warning('forceUpdate13DBFigureMetrics_R2FINAL: 输入 figure 无效。');
        return;
    end

    figure(fig);
    drawnow;

    axAll = findall(fig, 'Type', 'axes');

    % 排除 legend/colorbar 之类非普通 axes
    keep = true(numel(axAll),1);
    for i = 1:numel(axAll)
        tagi = string(get(axAll(i), 'Tag'));
        if contains(lower(tagi), 'legend') || contains(lower(tagi), 'colorbar')
            keep(i) = false;
        end
    end
    axAll = axAll(keep);

    pos = zeros(numel(axAll),4);
    for i = 1:numel(axAll)
        pos(i,:) = get(axAll(i), 'Position');
    end

    % 按从上到下、从左到右排序
    [~, ord] = sortrows([-pos(:,2), pos(:,1)]);
    axAll = axAll(ord);

    fprintf('\n================ R2 FINAL: 从当前 %s 图中散点强制重算 ================\n', modelName);

    dbCounter = 0;

    for ia = 1:numel(axAll)
        ax = axAll(ia);
        hScat = findall(ax, 'Type', 'Scatter');
        if isempty(hScat)
            continue;
        end

        dbCounter = dbCounter + 1;

        % 选择点数最多的 scatter，避免误选其它对象
        nPts = zeros(numel(hScat),1);
        for k = 1:numel(hScat)
            nPts(k) = numel(get(hScat(k), 'XData'));
        end
        [~, imax] = max(nPts);
        h = hScat(imax);

        x = get(h, 'XData');
        y = get(h, 'YData');
        x = x(:);
        y = y(:);

        valid = isfinite(x) & isfinite(y) & x > 1e-8 & ...
                x >= axisMin & x <= axisMax & y >= axisMin & y <= axisMax;
        x = x(valid);
        y = y(valid);

        N = numel(x);
        if N >= 2
            ARE = mean(abs((y - x) ./ x), 'omitnan');
            RMSE = sqrt(mean((y - x).^2, 'omitnan'));
            ssRes = sum((x - y).^2, 'omitnan');
            ssTot = sum((x - mean(x, 'omitnan')).^2, 'omitnan');
            if ssTot > eps
                R2 = 1 - ssRes / ssTot;
            else
                R2 = NaN;
            end
            Rtmp = corr(x, y, 'Rows', 'complete');
            R2corr = Rtmp.^2;
        else
            ARE = NaN; RMSE = NaN; R2 = NaN; R2corr = NaN;
        end

        fprintf('R2 FINAL | %s DB %02d: N=%d, ARE=%.5f, RMSE=%.5f, R2=%.5f, R2_corr=%.5f\n', ...
            modelName, dbCounter, N, ARE, RMSE, R2, R2corr);

        % 删除原有蓝色 R2 文本，避免两个 R2 重叠
        oldText = findall(ax, 'Type', 'Text');
        for it = 1:numel(oldText)
            sRaw = get(oldText(it), 'String');
            sText = string(sRaw);
            sHasR = any(contains(sText(:), 'R'));   % 标题是多行 cell/string 时，必须用 any 转成标量逻辑
            cText = get(oldText(it), 'Color');
            cText = cText(:);
            isBlueText = numel(cText) >= 3 && cText(3) > 0.5 && cText(1) < 0.3 && cText(2) < 0.3;
            if sHasR && isBlueText
                delete(oldText(it));
            end
        end

        title(ax, sprintf('DB %d, N=%d\nARE=%.3f, R^2=%.3f', dbCounter, N, ARE, R2), ...
            'FontName','Times New Roman', 'FontSize', 9.5, 'FontWeight','normal');

        text(ax, 0.05, 0.05, sprintf('R^2 = %.3f', R2), ...
            'Units','normalized', 'FontName','Times New Roman', ...
            'FontSize', 8.5, 'Color','b');
    end

    sgtitle(sprintf('%s: ML embedded drift-flux model in 13 databases  [R2 FINAL]', modelName), ...
        'FontName','Times New Roman', 'FontSize', 17, 'FontWeight','normal');

    drawnow;

    safeName = makeSafeFileNamePIEOneClick(modelName);
    exportgraphics(fig, fullfile(outFolder, sprintf('%s_13DB_Overview_R2FINAL.png', safeName)), 'Resolution', 300);
    savefig(fig, fullfile(outFolder, sprintf('%s_13DB_Overview_R2FINAL.fig', safeName)));

    fprintf('================ R2 FINAL: %s 图标题和蓝色 R2 已强制覆盖并重新保存 ================\n', modelName);
end

function safeName = makeSafeFileNamePIEOneClick(nameIn)

    safeName = char(string(nameIn));

    safeName = strrep(safeName, ' ', '_');
    safeName = strrep(safeName, '/', '_');
    safeName = strrep(safeName, '\', '_');
    safeName = strrep(safeName, ':', '_');
    safeName = strrep(safeName, '*', '_');
    safeName = strrep(safeName, '?', '_');
    safeName = strrep(safeName, '"', '_');
    safeName = strrep(safeName, '<', '_');
    safeName = strrep(safeName, '>', '_');
    safeName = strrep(safeName, '|', '_');
    safeName = strrep(safeName, '.', '');
    safeName = strrep(safeName, ',', '');

    while contains(safeName, '__')
        safeName = strrep(safeName, '__', '_');
    end
end
