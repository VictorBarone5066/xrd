from matplotlib import pyplot as plt

PLT_DATA_LOC = "plt.csv"

x, y = [], []
with open(PLT_DATA_LOC, 'r') as infile:
    for line in infile:
        lin = line.split(',')
        x.append(float(lin[0]))
        y.append(float(lin[1]))
    infile.close()

plt.plot(x, y, 'k', linewidth=0.8)
plt.xlim(10., 50.)
plt.ylim(0, 1.05)
plt.show()
