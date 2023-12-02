function output = u_fit(Lane_L_X, Lane_R_X, Lane_Y, inputIMG)
% 输入两列数据，分别进行二次函数拟合，将拟合结果以矩阵返回，并显示在输入的图像上

% x_L=a_L * x.^2 + b_L * x + c_L;
% x_R=a_R * x.^2 + b_R * x + c_R;

% 使用 polyfit 进行二次函数拟合
degree = 2; % 拟合的多项式阶数
coefficients_L = polyfit(Lane_Y, Lane_L_X, degree);
coefficients_R = polyfit(Lane_Y, Lane_R_X, degree);

% 从拟合系数中提取二次函数的系数
a_L = coefficients_L(1);
b_L = coefficients_L(2);
c_L = coefficients_L(3);
a_R = coefficients_R(1);
b_R = coefficients_R(2);
c_R = coefficients_R(3);

% 创建拟合的二次函数
Y=linspace(0,size(inputIMG,1));
X_L=a_L * Y.^2 + b_L * Y + c_L;
X_R=a_R * Y.^2 + b_R * Y + c_R;

lastlines=[399,400,401];
X_L_last=a_L * lastlines.^2 + b_L * lastlines + c_L;
X_R_last=a_R * lastlines.^2 + b_R * lastlines + c_R;
output=550-round((sum(X_R_last)+sum(X_L_last))/6);

%figure;
imshow(inputIMG);
hold on;
scatter(Lane_L_X, Lane_Y,'red','filled');
scatter(Lane_R_X, Lane_Y,'red','filled');
title('二次函数拟合');
plot(X_L, Y, 'r', 'LineWidth', 2, 'DisplayName', 'Left Lane'); % 绘制拟合曲线
plot(X_R, Y, 'r', 'LineWidth', 2, 'DisplayName', 'Right Lane');
scatter(550-output, 390,'green','filled');
hold off;

end