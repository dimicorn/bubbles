from functions import*


def main():
    a = str(input('Data or plot? '))
    b = ''
    c = ''
    if a == 'plot':
        b = str(input('Scatter or no? '))
        if b == 'no':
            c = str(input('Velocity or pressure? '))
    if a == 'data':
        data()
    if b == 'scat':
        scat()
    elif c == 'vel':
        vel_plots(2)
    elif c == 'pres':
        pres_plots(3)


if __name__ == '__main__':
    main()
