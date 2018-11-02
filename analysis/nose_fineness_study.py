import sys
sys.path.append('.')

from thermofluids import Atmosphere
from components import Rocket, Motor
from trajectory import TwoDTrajectorySolver
from aerodynamics import DragModel
from plotting import plotCurvesWithMax
import matplotlib.pyplot as plt
import numpy as np
import math
from sensitivity_analysis import sensitivity_analysis_apogee

atmosphere = Atmosphere(0)
launchAngle = 84*math.pi/180.0 #rad
launchrailLength = 8.0 #m

# Initial conditions
# y_0 is 0.01 instead of 0 so that the hit_ground event works. look into doing it another way
# x, vx, y, vy
y0 = [0., 0., 0.1, 0.]

rocketFilePath = "2019_rocket_fineness1.csv"
myRocket = Rocket(rocketFilePath)

motorFilePath = "Cesaroni_N5800.eng"
myMotor = Motor(motorFilePath)

solver = TwoDTrajectorySolver(
    myRocket, myMotor, atmosphere, y0, launchAngle, launchrailLength,\
    dragModel=DragModel.CoefficientsBox
)

solver.solve()
bApogee = solver.apogee

studyParameters = ["nose_length"]
start = 1
stop = 10
step = 1

counter = 1

apogeeSensitivityList = []
stabiltiySensitivityList = []
apogeeList = []
stabilityList = []

multipliers = np.arange(start, stop+step, step)

for multiplier in multipliers:
    print("({}/{}) {}".format(counter, len(multipliers), studyParameters[0]))
    counter += 1

    # start from baseline rocket every time
    rocketFilePath = "2019_rocket_fineness1.csv"
    myRocket = Rocket(rocketFilePath)

    motorFilePath = "Cesaroni_N5800.eng"
    myMotor = Motor(motorFilePath)
    
    myRocket.geometry["nose_length"] *= multiplier
    myRocket.geometry["body_length"] = 130/39.37 - myRocket.geometry["nose_length"] - myRocket.geometry["boattail_length"]
    myRocket.updateGeometry()
    
    solver = TwoDTrajectorySolver(
        myRocket, myMotor, atmosphere, y0, launchAngle, launchrailLength,\
        dragModel=DragModel.CoefficientsBox
    )
    
    solver.solve()
    
    apogeeList.append(solver.apogee)
    apogeeSensitivityList.append((solver.apogee - bApogee) / bApogee)


f, axarr = plt.subplots(2, sharex=True)

axarr[0].plot(multipliers, apogeeList, '-o')
axarr[0].set_title("apogee sensitivity with rocket parameters")
axarr[0].legend(studyParameters)

axarr[1].plot(multipliers, apogeeSensitivityList, '-o')
axarr[1].set_title("stability sensitivity with rocket parameters")
axarr[1].legend(studyParameters)

plt.show()