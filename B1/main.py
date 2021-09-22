from functions import*


def main():
    a = str(input('Data or plot? '))
    # b = ''
    c = ''
    d = ''
    if a == 'plot':
        b = str(input('Scatter or no? '))
        if b == 'scat':
            d = str(input('3d or 2d? '))
        elif b == 'no':
            c = str(input('Velocity or pressure? '))
    if a == 'data':
        data()
    if d == '3d':
        scat()
    elif d == '2d':
        colormap()
    elif c == 'vel':
        vel_plots(2)
    elif c == 'pres':
        pres_plots(3)


if __name__ == '__main__':
    main()
