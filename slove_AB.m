%%%%%%%%%%%%%%%%%%%%%%%%%Step 0：重置程序，定义变量%%%%%%%%%%%%%%%%%%%%%%%%%
tic
clear all
clc

% 定义机器人机体参数
syms R_w                % 驱动轮半径
syms R_l                % 驱动轮轮距/2
syms l_l l_r            % 左右腿长
syms l_wl l_wr          % 驱动轮质心到左右腿部质心距离
syms l_bl l_br          % 机体质心到左右腿部质心距离
syms l_c                % 机体质心到腿部关节中心点距离
syms m_w m_l m_b        % 驱动轮质量 腿部质量 机体质量
syms I_w                % 驱动轮转动惯量           (自然坐标系法向)
syms I_ll I_lr          % 驱动轮左右腿部转动惯量    (自然坐标系法向，实际上会变化)
syms I_b                % 机体转动惯量             (自然坐标系法向)
syms I_z                % 机器人z轴转动惯量        (简化为常量)

% 定义其他独立变量并补充其导数
syms theta_wl   theta_wr real % 左右驱动轮转角
syms dtheta_wl  dtheta_wr real
syms ddtheta_wl ddtheta_wr ddtheta_ll ddtheta_lr ddtheta_b real

% 定义状态向量
syms s ds phi dphi theta_ll dtheta_ll theta_lr dtheta_lr theta_b dtheta_b real

% 定义控制向量
syms T_wl T_wr T_bl T_br real

% 输入物理参数：重力加速度
syms g



%%%%%%%%%%%%%%%%%%%%%%%Step 1：解方程，求控制矩阵A，B%%%%%%%%%%%%%%%%%%%%%%%

% 通过原文方程组(3.11)-(3.15)，求出ddtheta_wl,ddtheta_wr,ddtheta_ll,ddtheta_lr,ddtheta_b表达式
eqn1 = (I_w*l_l/R_w+m_w*R_w*l_l+m_l*R_w*l_bl)*ddtheta_wl+(m_l*l_wl*l_bl-I_ll)*ddtheta_ll+(m_l*l_wl+m_b*l_l/2)*g*theta_ll+T_bl-T_wl*(1+l_l/R_w)==0;
eqn2 = (I_w*l_r/R_w+m_w*R_w*l_r+m_l*R_w*l_br)*ddtheta_wr+(m_l*l_wr*l_br-I_lr)*ddtheta_lr+(m_l*l_wr+m_b*l_r/2)*g*theta_lr+T_br-T_wr*(1+l_r/R_w)==0;
eqn3 = -(m_w*R_w*R_w+I_w+m_l*R_w*R_w+m_b*R_w*R_w/2)*ddtheta_wl-(m_w*R_w*R_w+I_w+m_l*R_w*R_w+m_b*R_w*R_w/2)*ddtheta_wr-(m_l*R_w*l_wl+m_b*R_w*l_l/2)*ddtheta_ll-(m_l*R_w*l_wr+m_b*R_w*l_r/2)*ddtheta_lr+T_wl+T_wr==0;
eqn4 = (m_w*R_w*l_c+I_w*l_c/R_w+m_l*R_w*l_c)*ddtheta_wl+(m_w*R_w*l_c+I_w*l_c/R_w+m_l*R_w*l_c)*ddtheta_wr+m_l*l_wl*l_c*ddtheta_ll+m_l*l_wr*l_c*ddtheta_lr-I_b*ddtheta_b+m_b*g*l_c*theta_b-(T_wl+T_wr)*l_c/R_w-(T_bl+T_br)==0;
eqn5 = ((I_z*R_w)/(2*R_l)+I_w*R_l/R_w)*ddtheta_wl-((I_z*R_w)/(2*R_l)+I_w*R_l/R_w)*ddtheta_wr+(I_z*l_l)/(2*R_l)*ddtheta_ll-(I_z*l_r)/(2*R_l)*ddtheta_lr-T_wl*R_l/R_w+T_wr*R_l/R_w==0;
[ddtheta_wl,ddtheta_wr,ddtheta_ll,ddtheta_lr,ddtheta_b] = solve(eqn1,eqn2,eqn3,eqn4,eqn5,ddtheta_wl,ddtheta_wr,ddtheta_ll,ddtheta_lr,ddtheta_b);

% 定义状态向量
syms s ds phi dphi theta_ll dtheta_ll theta_lr dtheta_lr theta_b dtheta_b real

% 定义控制向量
syms T_wl T_wr T_bl T_br real
% 通过计算雅可比矩阵的方法得出控制矩阵A，B所需要的全部偏导数
J_A = jacobian([ddtheta_wl,ddtheta_wr,ddtheta_ll,ddtheta_lr,ddtheta_b],[theta_ll,theta_lr,theta_b]);
J_B = jacobian([ddtheta_wl,ddtheta_wr,ddtheta_ll,ddtheta_lr,ddtheta_b],[T_wl,T_wr,T_bl,T_br]);

% 定义矩阵A，B，将指定位置的数值根据上述偏导数计算出来并填入
A = sym('A',[10 10])
B = sym('B',[10 4])

% 填入A数据：a25,a27,a29,a45,a47,a49,a65,a67,a69,a85,a87,a89,a105,a107,a109
for p = 5:2:9
    A_index = (p - 3)/2;
    A(2,p) = R_w*(J_A(1,A_index) + J_A(2,A_index))/2;
    A(4,p) = (R_w*(- J_A(1,A_index) + J_A(2,A_index)))/(2*R_l) - (l_l*J_A(3,A_index))/(2*R_l) + (l_r*J_A(4,A_index))/(2*R_l);
    for q = 6:2:10
        A(q,p) = J_A(q/2,A_index);
    end
end

% A的以下数值为1：a12,a34,a56,a78,a910，其余数值为0
for r = 1:10
    if rem(r,2) == 0
        A(r,1) = 0; A(r,2) = 0; A(r,3) = 0; A(r,4) = 0; A(r,6) = 0; A(r,8) = 0; A(r,10) = 0;
    else
        A(r,:) = zeros(1,10);
        A(r,r+1) = 1;
    end
end

% 填入B数据：b21,b22,b23,b24,b41,b42,b43,b44,b61,b62,b63,b64,b81,b82,b83,b84,b101,b102,b103,b104,
for h = 1:4
    B(2,h) = R_w*(J_B(1,h) + J_B(2,h))/2;
    B(4,h) = (R_w*(- J_B(1,h) + J_B(2,h)))/(2*R_l) - (l_l*J_B(3,h))/(2*R_l) + (l_r*J_B(4,h))/(2*R_l);
    for f = 6:2:10
        B(f,h) = J_B(f/2,h);
    end
end

% B的其余数值为0
for e = 1:2:9
    B(e,:) = zeros(1,4);
end



%%%%%%%%%%%%%%%%%%%%%Step 2：输入参数（可以修改的部分）%%%%%%%%%%%%%%%%%%%%%

% 物理参数赋值（唯一此处不可改变！），后面的数据通过增加后缀_ac区分模型符号和实际数据
g_ac = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      此处可以输入机器人机体基本参数                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%机器人机体与轮部参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % R1=0.0603;                         %驱动轮半径
    % L1=leg_length/2;                  %摆杆重心到驱动轮轴距离
    % LM1=leg_length/2;                 %摆杆重心到其转轴距离
    % l1=0.011;                          %机体质心距离转轴距离
    % mw1=0.6;                         %驱动轮质量
    % mp1=0.045;                         %杆质量
    % M1=1.44;                          %机体质量
    % Iw1=0.5*mw1*R1^2;                     %驱动轮转动惯量
    % Ip1=mp1*((L1+LM1)^2+0.048^2)/12.0; %摆杆转动惯量
    % IM1=M1*(0.135^2+0.066^2)/12.0;       %机体绕质心转动惯量

% mp1=0.045;                                 %杆质量
mp1=0.06; 
R_w_ac = 0.05;                           % 驱动轮半径                  （单位：m）
R_l_ac = 0.155;                            % 两个驱动轮之间距离/2         （单位：m）
l_c_ac = 0.00001;                            % 机体质心到腿部关节中心点距离  （单位：m）
m_w_ac = 0.4;                              % 驱动轮质量
% m_l_ac = 0.124;    %什么东西            % 腿部质量 
m_l_ac = 0.06;
m_b_ac = 2;                             % 机体质量  （单位：kg）
I_w_ac = 0.5*m_w_ac*R_w_ac^2;                        % 驱动轮转动惯量               （单位：kg m^2）
I_b_ac = m_b_ac*(0.3^2+0.3^2)/12.0;                  % 机体转动惯量(自然坐标系法向)  （单位：kg m^2）
% I_z_ac = 0.226;       %什么东西                     % 机器人z轴转动惯量            （单位：kg m^2）
I_z_ac = 0.0314415;
%%%%%%%%%%%%%%%%%%%%%%机器人腿部参数（定腿长调试用）%%%%%%%%%%%%%%%%%%%%%%%%

% 如果需要使用此部分，请去掉120-127行以及215-218行注释，然后将224行之后的所有代码注释掉
% 或者点击左侧数字"224"让程序只运行行之前的内容并停止

leg_length = 0.218;
l_l_ac = 0.218;        % 左腿摆杆长度                      （左腿对应数据）  （单位：m）
l_wl_ac = leg_length/2;       % 左驱动轮质心到左腿摆杆质心距离                      （单位：m）
l_bl_ac = leg_length/2;       % 机体转轴到左腿摆杆质心距离                          （单位：m）
% I_ll_ac = 0.02054500;    %什么东西  % 左腿摆杆转动惯量                                   （单位：kg m^2）
l_r_ac = 0.218;        % 右腿摆杆长度                      （右腿对应数据）  （单位：m）
l_wr_ac = leg_length/2;       % 右驱动轮质心到右腿摆杆质心距离                      （单位：m）
l_br_ac = leg_length/2;       % 机体转轴到右腿摆杆质心距离                          （单位：m）
% I_lr_ac = 0.02054500;      % 右腿摆杆转动惯量                                   （单位：kg m^2）
I_lr_ac=mp1*((l_wr_ac+l_br_ac)^2+0.01^2)/12.0;   %右腿摆杆转动惯量
I_ll_ac=mp1*((l_wl_ac+l_bl_ac)^2+0.01^2)/12.0;   %左腿摆杆转动惯量

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 3：代入实际参数并数值化 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 若矩阵中仍包含角度变量（theta_ll、theta_lr、theta_b），
% 先在平衡点（小扰动线性化点）代入零姿态
A_lin = subs(A, [theta_ll, theta_lr, theta_b], [0, 0, 0]);
B_lin = subs(B, [theta_ll, theta_lr, theta_b], [0, 0, 0]);

% 将所有物理参数代入实际值
A_subs = subs(A_lin, ...
    [R_w, R_l, l_l, l_r, l_wl, l_wr, l_bl, l_br, l_c, ...
     m_w, m_l, m_b, I_w, I_ll, I_lr, I_b, I_z, g], ...
    [R_w_ac, R_l_ac, l_l_ac, l_r_ac, l_wl_ac, l_wr_ac, l_bl_ac, l_br_ac, l_c_ac, ...
     m_w_ac, m_l_ac, m_b_ac, I_w_ac, I_ll_ac, I_lr_ac, I_b_ac, I_z_ac, g_ac]);

B_subs = subs(B_lin, ...
    [R_w, R_l, l_l, l_r, l_wl, l_wr, l_bl, l_br, l_c, ...
     m_w, m_l, m_b, I_w, I_ll, I_lr, I_b, I_z, g], ...
    [R_w_ac, R_l_ac, l_l_ac, l_r_ac, l_wl_ac, l_wr_ac, l_bl_ac, l_br_ac, l_c_ac, ...
     m_w_ac, m_l_ac, m_b_ac, I_w_ac, I_ll_ac, I_lr_ac, I_b_ac, I_z_ac, g_ac]);

% 将符号矩阵转换为数值矩阵（保留6位小数）
A_num = double(vpa(A_subs, 6));
B_num = double(vpa(B_subs, 6));

% 输出结果
disp('================ A矩阵（数值化后） ================');
disp(A_num);

disp('================ B矩阵（数值化后） ================');
disp(B_num);

% 保存结果以便后续导入Python或控制器设计
save('system_matrices.mat', 'A_num', 'B_num');

toc
