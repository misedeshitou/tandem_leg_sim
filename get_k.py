import sympy as sp
import numpy as np
from scipy.linalg import solve_continuous_are
from control import lqr
from scipy.io import loadmat
from scipy.io import savemat

# ========================= Step 1: 定义符号变量 =========================
R_w, R_l = sp.symbols('R_w R_l', real=True)
l_l, l_r, l_wl, l_wr, l_bl, l_br, l_c = sp.symbols('l_l l_r l_wl l_wr l_bl l_br l_c', real=True)
m_w, m_l, m_b = sp.symbols('m_w m_l m_b', real=True)
I_w, I_ll, I_lr, I_b, I_z = sp.symbols('I_w I_ll I_lr I_b I_z', real=True)
g = sp.Symbol('g', real=True)

theta_ll, theta_lr, theta_b = sp.symbols('theta_ll theta_lr theta_b', real=True)
T_wl, T_wr, T_bl, T_br = sp.symbols('T_wl T_wr T_bl T_br', real=True)
ddtheta_wl, ddtheta_wr, ddtheta_ll, ddtheta_lr, ddtheta_b = sp.symbols('ddtheta_wl ddtheta_wr ddtheta_ll ddtheta_lr ddtheta_b', real=True)

# 读取 .mat 文件
data = loadmat('system_matrices.mat')
# data= loadmat('AB.mat')

# 提取 A, B
A = data['A_num']
B = data['B_num']
# A = data['A']
# B = data['B']

print("A =", A)
print("B =", B)

# ========================= Step 6: 定义 Q, R =========================
#上交建模
#        s     ds     phi     dphi     theta_ll dtheta_ll theta_lr dtheta_lr theta_b dtheta_b
Q = np.array([
    [1,   0,   0,  0,   0,   0,   0,  0,   0,   0],
    [0,   5,   0,  0,   0,   0,   0,  0,   0,   0],
    [0,   0,   1,  0,   0,   0,   0,  0,   0,   0],
    [0,   0,   0,  5,   0,   0,   0,  0,   0,   0],
    [0,   0,   0,  0, 150,   0,   0,  0,   0,   0],
    [0,   0,   0,  0,   0,   1,   0,  0,   0,   0],
    [0,   0,   0,  0,   0,   0, 150,  0,   0,   0],
    [0,   0,   0,  0,   0,   0,   0,  1,   0,   0],
    [0,   0,   0,  0,   0,   0,   0,  0, 300,   0],
    [0,   0,   0,  0,   0,   0,   0,  0,   0,  10]
], dtype=float)

R = np.array([
    [10, 0,  0,  0],
    [ 0,10,  0,  0],
    [ 0, 0, 20,  0],
    [ 0, 0,  0, 20]
], dtype=float)


# 哈工程建模
# Q = np.array([
#     [150, 0, 0, 0, 0, 0],
#     [0, 1, 0, 0, 0, 0],
#     [0, 0, 1, 0, 0, 0],
#     [0, 0, 0, 10, 0, 0],
#     [0, 0, 0, 0, 300, 0],
#     [0, 0, 0, 0, 0, 10]
# ], dtype=float)

# R = np.array([
#     [90, 0],
#     [0, 5]
# ], dtype=float)
# ========================= Step 7: LQR 求解 =========================

K, S, E = lqr(A, B, Q, R)

print("LQR 增益矩阵 K：")
print(np.array2string(K, formatter={'float_kind':lambda x: f"{x: .5f}"}))

savemat('system_matrices_from_py.mat', {'A': A, 'B': B, 'K': K})