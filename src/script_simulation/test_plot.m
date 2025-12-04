H = diag([-2,-2,-2]);
f = [
    5.9628
   22.7806
   13.9344
   ];
c = -91.6592;
range_val = 100;

plot_quadratic(H,f,c,range_val);


function plot_quadratic(H, f, c, range_val)
    % plot_quadratic: 이차 형식 x'Hx + f'x + c를 플롯합니다.
    % H: 2x2 행렬
    % f: 2x1 또는 1x2 벡터
    % c: 스칼라
    % range_val: 플롯할 x1, x2의 범위 (예: 10이면 -10~10)

    if nargin < 4
        range_val = 10; % 기본 범위 설정
    end

    % 1. 그리드 생성 (Meshgrid)
    step = range_val / 50; % 해상도 조절
    [X1, X2] = meshgrid(-range_val:step:range_val, -range_val:step:range_val);
    
    % 2. 함수 값 계산 (벡터화 연산)
    % 수식: Z = H(1,1)x1^2 + (H(1,2)+H(2,1))x1x2 + H(2,2)x2^2 + f1x1 + f2x2 + c
    
    term_quad = H(1,1).*X1.^2 + (H(1,2) + H(2,1)).*X1.*X2 + H(2,2).*X2.^2;
    term_lin  = f(1).*X1 + f(2).*X2;
    Z = term_quad + term_lin + c;

    % 3. 시각화 (서브플롯 활용)
    figure('Color', 'w', 'Position', [100, 100, 1000, 500]);
    
    % [왼쪽] 3차원 서피스 플롯
    subplot(1, 2, 1);
    surfc(X1, X2, Z); % surfc는 바닥에 등고선을 같이 그려줍니다
    shading interp;   % 색상을 부드럽게
    colormap jet;
    colorbar;
    title('3D Surface Plot ($x^T H x + f^T x + c$)', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('x_1'); ylabel('x_2'); zlabel('f(x)');
    grid on;
    axis square;

    % [오른쪽] 2차원 등고선 플롯 (위에서 본 모습)
    subplot(1, 2, 2);
    contourf(X1, X2, Z, 20); % 20개의 등고선 레벨
    colorbar;
    title('2D Contour Plot', 'FontSize', 14);
    xlabel('x_1'); ylabel('x_2');
    grid on;
    axis square;
    
    % 4. 극점(Stationary Point) 표시 (H가 역행렬을 가질 때)
    if det(H) ~= 0
        x_opt = -H \ f; % x* = -H^-1 * f
        % 범위 안에 있을 때만 표시
        if abs(x_opt(1)) <= range_val && abs(x_opt(2)) <= range_val
            z_opt = x_opt' * H * x_opt + f' * x_opt + c;
            
            subplot(1, 2, 1); hold on;
            plot3(x_opt(1), x_opt(2), z_opt, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
            text(x_opt(1), x_opt(2), z_opt, '  Optimum', 'Color', 'r', 'FontWeight', 'bold');
            
            subplot(1, 2, 2); hold on;
            plot(x_opt(1), x_opt(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
            text(x_opt(1), x_opt(2), '  Optimum', 'Color', 'r', 'FontWeight', 'bold');
        end
    end
end