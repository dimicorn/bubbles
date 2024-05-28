# All values for gamma
G = [1.01, 1.05, 1.1, 1.4]
# Same for k_rho
K = [0, 1, 1.5, 2, 2.5, 3]
# And for n_int too
N_int = [0.5, 1, 1.5, 2.5]

# Astronomical unit in meters
au = 150 * 10**9

# seconds in 1 year
t_y = 3600 * 24 * 365

# Starting and values for
# velocity
v0 = 1000
vk = 100000

# number of years
numb0 = 1000
numb_k = 100000

# ratio of R_c/R_sw
k0 = 10
k_k = 1000

# Place where the Excel table wih all the data will be created
file_path = "C:/Users/dmitr/Desktop/MIPT/Лето 21/Dr. Shang/Bubbles/Data/Auto2.xlsx"

if __name__ == '__main__':
    print('You are running the wrong file!')
