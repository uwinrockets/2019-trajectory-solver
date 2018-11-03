n = 4

for i in range(0, n):
    for j in range(1, n-i):
        print("{} {}".format(i, j))

from pyDOE import *
np.savetxt('name.csv', lhs(4, samples=100))