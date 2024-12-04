clc;clear;
rng(47);TN = 10;T1 = 1;t0 = 0.1;j = 13;
delta = 2.^[-j+6, -j+5, -j+4, -j+3, -j];
N = round(TN / delta(end));    % 模拟多个样本数量
Nu = round(T1 / delta(end)) + 1;  % 供使用的样本数量
t = t0:delta(end):TN;    % 时间轴
alpha_list = 0.3:0.1:0.9;
num_iterations = 1000;  % 样本数量
kappa = 8;sigma = 0.5;theta = 0.125;
theta_adjusted = theta - sigma^2/(4*kappa);  
convergence_rate = zeros(2, length(alpha_list));
convergence_rate(1, :) = alpha_list;
delta_errors = zeros(length(delta) - 1, num_iterations);  % 误差矩阵
x_solutions = zeros(length(delta), num_iterations);       % 解矩阵

for alpha_idx = 1:length(alpha_list)    % 遍历alpha
	alpha = alpha_list(alpha_idx);
	phase_shift = pi / 2;               
	
	for path_idx = 1:num_iterations    % 遍历路径
		% 模拟稳定过程 D
		stable_process = zeros(1, N + 1);    % 稳定过程矩阵
		for i = 1:N
			random_angle = (rand - 0.5) * pi;     
			exp_random = exprnd(1);               
			stable_increment = sin(alpha * (random_angle + phase_shift)) / ...
							cos(random_angle)^(1/alpha) * ...
							(cos(random_angle - alpha * (random_angle + phase_shift)) / exp_random)^((1 - alpha) / alpha);
			stable_process(i + 1) = stable_process(i) + delta(end)^(1/alpha) * stable_increment;
		end
		D = stable_process;

		% 模拟布朗运动
		brownian_increments = sqrt(delta(end)) * randn(1, N);    
		brownian_motion = cumsum(brownian_increments);           

		% 找到D中第一次超过1的位置
		index = find(D > 1, 1, 'first');
		D = D(1:index);

		% 模拟E
		E = zeros(1, Nu);
		tt = 0:delta(end):T1;
		for i = 1:length(tt)
			E(i) = find(D > tt(i), 1, 'first') * delta(end);
		end

		E_list = zeros(length(delta), Nu);
		for i = 1:length(delta)
			step = round(delta(i) / delta(end));
			E_list(i, 1:round(T1 / delta(i)) + 1) = E(1:step:end);
		end


		for z = 1:length(delta)
			X = ones(1,T1/delta(z)+1);
			for i = 1:round(T1 / delta(z))
				Winc = brownian_motion(E_list(z, i + 1) / delta(end)) - brownian_motion(E_list(z, i) / delta(end));
				Einc = E_list(z, i + 1) - E_list(z, i);

				if Einc == 0
					X(i + 1) = X(i);
				else
					err = 1;
					while err > 1e-6
						Xtemp=X(i+1);
						f = X(i+1) - X(i) -0.5*kappa*(theta_adjusted/X(i+1) - X(i+1))*Einc - 0.5*sigma*Winc;
						f1 = 1 + 0.5*kappa*Einc + kappa*theta_adjusted*Einc/(2*X(i+1)^2);
						X(i + 1) = X(i+1) - f / f1;
						err = abs(X(i + 1) - Xtemp);
					end
				end
			end
			Xbem(z) = X(end);
		end
		x_solutions(:,path_idx) = Xbem;
		fprintf('alpha = %d, 正在计算第 %d 条路径\n', alpha_list(alpha_idx), path_idx);
	end

	% 计算误差和收敛率
	for delta_idx = 1:length(delta) - 1
		delta_errors(delta_idx, :) = abs(x_solutions(delta_idx, :) - x_solutions(end, :));
	end
	mean_errors = mean(delta_errors, 2)';
	delta_values = delta(1:end-1);
	[correlation, slope, intercept] = regression(log(delta_values), log(mean_errors));
	convergence_rate(2, alpha_idx) = slope;
end
convergence_rate

