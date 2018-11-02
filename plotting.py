import matplotlib.pyplot as plt
import numpy as np

# change so it returns a plot instead
def plotCurvesWithMax(x, y, labelStr, title):
  plt.plot(x, y, label=labelStr)
  plt.title(title)
  
  apogee_y = np.max(y)
  apogee_x= x[np.argmax(y)]
  
  plt.text(apogee_x, apogee_y+(0.02*apogee_y), "({:.2f}, {:.2f})".format(apogee_x, apogee_y))
  x1,x2,y1,y2 = plt.axis()
  plt.axis((x1,x2,y1,y2 + (0.05*apogee_y)))
  plt.legend()
  plt.show()