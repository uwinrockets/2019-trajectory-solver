from trajectory import TwoDTrajectorySolver
from components import Rocket, Motor
import pprint

def sensitivity_analysis_apogee(parameterKey, solver, multipliers):
    
    solver.solve()
    baselineApogee = solver.apogee
    baselineStability = solver.rocket.staticStability

    baselineParameter = float(solver.rocket.geometry[parameterKey])
    apogeeSensitivityList = []
    newApogeeList = []
    stabilitySensitivityList = []
    newStabilityList = []

    counter = 0

    for multiplier in multipliers:
        newParameter = baselineParameter * (1 + multiplier)
        solver.rocket.geometry[parameterKey] = newParameter
        solver.rocket.updateGeometry()

        solver.solve()
        
        newApogee = solver.apogee
        apogee_sensitivity = ((newApogee - baselineApogee) / baselineApogee)
        apogeeSensitivityList.append(apogee_sensitivity)
      
        newStability = solver.rocket.staticStability
        stability_sensitivity = ((newStability - baselineStability) / baselineStability)
        stabilitySensitivityList.append(stability_sensitivity)

        newApogeeList.append(newApogee)
        newStabilityList.append(newStability)

        counter += 1
        print("Iteration: {}/{}".format(counter, len(multipliers)))

    results = {
        'baseline_parameter': baselineParameter,
        'baseline_apogee': baselineApogee,
        'baseline_sensitivity': baselineStability,
        'apogee_sensitivity': apogeeSensitivityList,
        'stability_sensitivity': stabilitySensitivityList,
        'apogee_list': newApogeeList,
        'stability_list': newStabilityList,
    }

    return results