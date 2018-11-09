from __future__ import division
from scipy import interpolate
import numpy as np
import csv
import math

class Motor:
  '''
  Motor class for a single stage rocket
  Assumes that time=0 is ignition (largely simplifies things)
  '''
  def __init__(self, filepath):
    self.filepath = filepath
    self.motor_data = np.loadtxt(self.filepath, skiprows=2)
    self.time_data = self.motor_data[:,0]
    self.thrust_data = self.motor_data[:,1]
    self.thrust_fit = interpolate.interp1d(self.time_data, self.thrust_data, kind="linear", bounds_error=False, fill_value=0)
    self.wetMass = None
    self.dryMass = None
    self.burnTime = self.time_data[-1]
 
    with open(self.filepath, 'r') as f:
      first_line = f.readline().split()
      self.dryMass = float(first_line[1])
      self.wetMass = float(first_line[2])
      self.massFlowRate = (self.wetMass - self.dryMass) / self.burnTime
      
  def getThrust(self, time):
    return self.thrust_fit(time)
  
  def isPowered(self, time):
    if (time >= self.burnTime) or (time<0):
      return False
    else:
      return True
  
  def getMass(self, time):
    if time < 0:
      return self.wetMass
    elif time > self.burnTime:
      return self.dryMass
    else:
      return self.wetMass - self.massFlowRate * time
  
  def getMassFlow(self, time):
    if time < 0:
      return 0.0
    elif time > self.burnTime:
      return 0.0
    else:
      return self.massFlowRate

class Rocket:
    def __init__(self, rocketFilePath):
        with open(rocketFilePath) as rocketFile:
            reader = csv.reader(rocketFile)
            self.geometry = dict(reader)

            self.length_units = self.geometry["unit_system_length"]
            self.mass_units = self.geometry["unit_system_mass"]
            self.angle_units = "degrees"

            self.updateGeometry()

    # this is sketchy but i didnt think far ahead when i coded this in the first place
    # fix later to only reference the keys
    def updateGeometry(self):

      self.conversionFactorLength = 1.
      self.conversionFactorMass = 1.
      self.conversionFactorAngle = 1.

      if self.length_units == "m": #meters
          self.conversionFactorLength = 1.
      elif self.length_units == "mm": #mm
          self.conversionFactorLength = 1./1000.
      elif self.length_units == "ft": #feet
          self.conversionFactorLength = 1./3.281
      elif self.length_units == "in": #inches
          self.conversionFactorLength = 1./39.37

      if self.mass_units == "kg":
          self.conversionFactorMass = 1.
      elif self.mass_units == "g":
          self.conversionFactorMass = 1000.
      elif self.conversionFactorMass == "lbs":
          self.conversionFactorMass = 1./2.205

      if self.angle_units == "degrees":
          self.conversionFactorAngle = math.pi / 180.
      elif self.angle_units == "radians":
          self.conversionFactorAngle = 1

      self.lengthKeys = [
        "nose_diameter", "nose_length",
        "body_diameter", "body_length",
        "boattail_length", "boattail_diameter",
        "motor_diameter", "finThickness",
        "finRootChord", "finTipChord",
        "finSpan", "finSweepLength", "finConnection"
      ]
      
      self.massKeys = [
        "motorless_mass"
      ]

      self.angleKeys = [
        "finLESweep", "finLETotal"
      ]

      # This will convert everything into SI units
      for key in self.lengthKeys:
        if key in self.geometry:
          self.geometry[key] = float(self.geometry[key]) * self.conversionFactorLength
      self.length_units = "m"
      
      for key in self.massKeys:
        if key in self.geometry:
          self.geometry[key] = float(self.geometry[key]) * self.conversionFactorMass
      self.mass_units = "kg"
      
      for key in self.angleKeys:
        if key in self.geometry:
          self.geometry[key] = float(self.geometry[key]) * self.conversionFactorAngle
      self.angle_units = "radians"

      # single body diameter. remove in future versions to account for different sized diameters
      self.geometry["nose_diameter"] = self.geometry["body_diameter"]
      
      self.geometry["nFins"] = float(self.geometry["nFins"])

      self.geometry["referenceArea"] = math.pi / 4. * math.pow(self.geometry["body_diameter"], 2.)
      
      # self.geometry["taper"] = self.geometry["finTipChord"] / self.geometry["finRootChord"]
      
      #self.geometry["meanAeroChord"] = (2./3.) * self.geometry["finRootChord"] \
      #  * (1 + self.geometry["taper"] + self.geometry["taper"]**2) / (1 + self.geometry["taper"])
      
      #self.geometry["finAspectRatio"] = 2*self.geometry["finSpan"]\
      #  / ((1 + self.geometry["taper"]) * self.geometry["finRootChord"])
      
      #self.geometry["finArea"] = self.geometry["finSpan"]**2 / self.geometry["finAspectRatio"]
    
      self.geometry["totalLength"] = self.geometry["nose_length"]\
        + self.geometry["boattail_length"] + self.geometry["body_length"]

      if "finConnection" not in self.geometry:
        self.geometry["finConnection"] = self.geometry["nose_length"] + self.geometry["body_length"] - self.geometry["finRootChord"]
    
      if self.geometry["finSpan"] == 0:
        self.geometry["finLESweep"] = 0
      else:
        self.geometry["finLESweep"] = math.atan(self.geometry["finSweepLength"]/self.geometry["finSpan"])

      self.geometry["finMidChord"] = math.sqrt(self.geometry["finSpan"]**2\
        + math.pow(self.geometry["finSweepLength"] + self.geometry["finTipChord"]/2 - self.geometry["finRootChord"]/2, 2))

      self.centerGravity = 0.625 * self.geometry["totalLength"] # just for now. i got this from openrocket

      self.calculateStaticStability()
    
    def calculateStaticStability(self):
      # c_na_b = 2 # per radian
      # x_ac_b = 0.63 * self.geometry["nose_length"]
      
      # c_na_t = math.pi * self.geometry["finAspectRatio"] / 2
      # x_ac_t = self.geometry["finConnection"] + 0.25 * self.geometry["meanAeroChord"]
      
      # c1 = c_na_b * (self.centerGravity - x_ac_b) / self.geometry["body_diameter"]
      # c2 = c_na_t * (self.centerGravity - x_ac_t) / self.geometry["body_diameter"]
      
      # self.staticStability = -1 * (c1 + c2*(self.geometry["finArea"]/self.geometry["referenceArea"]))\
      # / (c_na_b + c_na_t * (self.geometry["finArea"]/self.geometry["referenceArea"]))

      c_t = self.geometry["finTipChord"]
      c_r = self.geometry["finRootChord"]
      l_m = self.geometry["finMidChord"]
      s = self.geometry["finSpan"]
      d = self.geometry["body_diameter"]
      r = d/2
      d_n = self.geometry["nose_diameter"]
      d_b = self.geometry["boattail_diameter"]
      n = self.geometry["nFins"]
      x_f = self.geometry["finConnection"]

      # nose
      c_na_n = 2
      x_ac_n = 0.466 * self.geometry["nose_length"]

      # boattail
      if (self.geometry["boattail_length"] == 0) or (self.geometry["boattail_diameter"] == 0):
        c_na_c = 0
        x_ac_c = 0
      else:
          c_na_c = 2 * (math.pow(d_b/d_n, 2) - math.pow(d/d_n, 2))
          x_ac_c = self.geometry["nose_length"] + self.geometry["body_length"]\
            * (self.geometry["boattail_length"]/3)*(1+(1-d/d_b)/(1-math.pow(d/d_b, 2)))

      # fins
      if (l_m == 0) or (c_r == 0) or (s == 0):
        c_na_f = 0
        x_ac_f = 0
      else:
        c_na_f = (1 + r/(s+r)) * (4*n*math.pow(s/d, 2))/(1+math.sqrt(1+math.pow((2*l_m)/(c_r+c_t), 2)))
        x_ac_f = x_f + (l_m/3)*(c_r+2*c_t)/(c_r+c_t) + (1/6)*((c_r+c_t) - (c_r*c_t)/(c_r+c_t))

      c_na_r = c_na_f + c_na_n + c_na_c

      self.centerPressure = (c_na_n*x_ac_n + c_na_f*x_ac_f + c_na_c*x_ac_c)/c_na_r
      self.staticStability = (self.centerPressure - self.centerGravity) / self.geometry["body_diameter"]
      
      # "error" checking
      for key in self.lengthKeys:
        if (self.geometry[key] < 0):
          self.staticStability = np.NINF