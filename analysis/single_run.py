from __future__ import division

import sys
sys.path.append('.')

from thermofluids import Atmosphere
from components import Rocket, Motor
from trajectory import TwoDTrajectorySolver
from aerodynamics import DragModel
from plotting import plotCurvesWithMax
from sensitivity_analysis import sensitivity_analysis_apogee

import matplotlib.pyplot as plt
import numpy as np
import math

# additional parameters the solver needs
atmosphere = Atmosphere(0)
launchAngle = 84*math.pi/180.0 #rad
launchrailLength = 8.0 #m

# Initial conditions
# y_0 is 0.01 instead of 0 so that the hit_ground event works. look into doing it another way
# x, vx, y, vy
y0 = [0., 0., 0.1, 0.]

rocketFilePath = "2019_rocket_genetic.csv"
myRocket = Rocket(rocketFilePath)

motorFilePath = "Cesaroni_N5800.eng"
myMotor = Motor(motorFilePath)

solver = TwoDTrajectorySolver(
    myRocket, myMotor, atmosphere, y0, launchAngle, launchrailLength,\
    dragModel=DragModel.CoefficientsBox
)

solver.solve()
solver.post_process()
print("Static Stability: {:.2f}".format(solver.rocket.staticStability))
solver.printPostProcess()

plt.plot(solver.t, solver.y, label="apogee")
plt.legend()
plt.show()

plt.plot(solver.t, solver.machVsTime, label="mach")
for key in solver.dragCoefficientsVsTime:
    plt.plot(solver.t, solver.dragCoefficientsVsTime[key], label=key)
plt.xlim(0, solver.apogee_time + 5)
plt.ylim(0, solver.maxMach + 0.1)
plt.legend()
plt.show()