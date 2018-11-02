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

rocketFilePath = "2019_rocket.csv"
myRocket = Rocket(rocketFilePath)

motorFilePath = "Cesaroni_N5800.eng"
myMotor = Motor(motorFilePath)

solver = TwoDTrajectorySolver(
    myRocket, myMotor, atmosphere, y0, launchAngle, launchrailLength,\
    dragModel=DragModel.CoefficientsBox
)


parameterKeys = ['finSpan', 'finTipChord', 'finRootChord', 'finSweepLength']
baselineParameters = {}
doe_coeff = np.loadtxt('analysis/lhc_n4_200samples_random.csv')
apogeeList = []
stabilityList = []

for key in parameterKeys:
    baselineParameters[key] = myRocket.geometry[key]

counter = 0

for row in doe_coeff:
    # change from baseline
    for i in range(len(parameterKeys)):
        myRocket.geometry[parameterKeys[i]] = baselineParameters[parameterKeys[i]] * row[i]
    myRocket.updateGeometry()
    
    solver.solve()
    apogeeList.append(solver.apogee)
    stabilityList.append(myRocket.staticStability)

    # set back to baseline
    for i in range(len(parameterKeys)):
        myRocket.geometry[parameterKeys[i]] = baselineParameters[parameterKeys[i]]

    counter += 1
    print("... Iteration {}/{}. Apogee {:.0f}m Stability {:.2f}".format(counter, len(doe_coeff), solver.apogee, myRocket.staticStability))

np.savetxt('analysis/apogee_list.csv', apogeeList)
np.savetxt('analysis/stability_list.csv', stabilityList)