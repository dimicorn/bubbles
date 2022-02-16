import matplotlib.pyplot as plt


x = []
y = []
file = open("data.txt", "r")
for i in file:
    print(i)
    j = i.split(" ")
    print(j[0], j[1])
    x.append(float(j[0]))
    y.append(float(j[1]))

plt.plot(x, y)
plt.savefig("data.png")
