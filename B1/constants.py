# All values of gamma
G = [1.01, 1.05, 1.1, 1.4, 5/3]

# Same for k_rho
K = [0, 1, 1.5, 2, 2.5]  # 3

# And n_int too
N_int = [0.5, 1, 1.5, 2.5]

# Colors for plot option
Colors = ['purple', 'blue', 'red', 'green', 'cyan', 'orange']

# Starting and ending values of lambda
x0 = 1
x_k = 0.8

# Initial conditions
v0 = 1
rho0 = 1
p0 = 1

# Plot labels
Gamma = 'gamma = '
K_rho = 'k_rho = '
N = 'n_int = '

# The place where will be created an excel table with all the data
file_path = "C:/Users/dmitr/Desktop/MIPT/Лето 21/Dr. Shang/Bubbles/Data/Auto1.xlsx"

au = 150 * 10**9

# For data option starting values of the parameters are being used
# Starting and ending values of velocity
v_int = 1000
vk = 100000

# Starting and ending number of years for time
n0 = 1000
nk = 100000

# Starting and ending values of ratio R_c/R_sw
k0 = 10
k_k = 1000

if __name__ == '__main__':
    print('You are running the wrong file!')
