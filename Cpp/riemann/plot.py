import sys
import matplotlib.pyplot as plt


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

def main():
    if sys.argv:
        print('No input file was given!')
        exit(0)
    
    input_file = sys.argv[1]
    coord, density, velocity, pressure, int_energy = parse(input_file)
    values = [density, velocity, pressure, int_energy]
    names = ['den', 'vel', 'pres', 'ener']
    labels = ['\\rho', 'v', 'p', 'U']

    for value, name, label in zip(values, names, labels):
        fig, ax = plt.subplots()
        ax.scatter(coord, value, marker='.')
        ax.set_xlabel('x')
        ax.set_ylabel(rf'${label}$')
        ax.set_title(rf'${label}(x)$')

        plt.savefig(f'figures/{name}.png')
        plt.close(fig)

if __name__ == '__main__':
    main()