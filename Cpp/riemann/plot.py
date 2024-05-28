import sys
import matplotlib.pyplot as plt
import os
import numpy as np


def parse(file_name: str) -> tuple[list, list, list, list, list]:
    coord, density, velocity = [], [], []
    pressure, int_energy = [], []
    with open(file_name) as f:
        contents = f.readlines()
        for line in contents:
            line = line.split()
            x, den, vel, pres, energy = tuple(map(float, line))
            coord.append(x)
            density.append(den)
            velocity.append(vel)
            pressure.append(pres)
            int_energy.append(energy)

    return (coord, density, velocity, pressure, int_energy)

def main() -> None:
    if len(sys.argv) == 1:
        print('No input file was given!')
        exit(1)
    
    input_file, gamma, dl, vl, pl = tuple(sys.argv[1:6])
    coord, density, velocity, pressure, int_energy = parse(input_file)
    # sound = np.array(np.sqrt(gamma * np.array(pressure)))
    values = [density, velocity, pressure, int_energy]
    names = ['den', 'vel', 'pres', 'ener'] #  'sound'
    labels = ['\\rho', 'v', 'p', 'U'] # 'a'
    if not os.path.exists(f'shang/gamma_{gamma}_dl_{dl}_pl_{pl}'):
            os.mkdir(f'shang/gamma_{gamma}_dl_{dl}_pl_{pl}')
    for value, name, label in zip(values, names, labels):
        fig, ax = plt.subplots()
        ax.scatter(coord, value, marker='.')
        ax.set_xlabel('x')
        if name == 'den':
            ax.set_yscale('log')
            ax.set_ylabel(rf'$\log({label})$')
        else:
            ax.set_ylabel(rf'${label}$')
        ax.set_title(rf'${label}(x), \gamma = {gamma}, \rho_L = {dl}, v_L = {vl}, p_L = {pl}$', loc='center')
        
        plt.savefig(f'shang/gamma_{gamma}_dl_{dl}_pl_{pl}/{name}.png')
        plt.close(fig)
    return 
    value, name, label = density, 'den', '\\log(\\rho)'
    fig, ax = plt.subplots()
    ax.scatter(coord, value, marker='.')
    ax.set_xlabel('x')
    ax.set_ylabel(rf'${label}$')
    ax.set_yscale('log')
    # ax.set_title(rf'${label}(x)$', loc='center')
    ax.set_title(rf'$\gamma = {gamma}, \rho_L = {dl}, p_L = {pl}$', loc='center')
    plt.savefig(f'shang_2/log_gamma_{gamma}_dl_{dl}_pl_{pl}.png')
    plt.close(fig)

if __name__ == '__main__':
    main()