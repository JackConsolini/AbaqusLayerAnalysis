""" 

Abaqus-python code for extracting state variable or layer thickness information from 2D multi-
layered simulations. Analysis is conducted based on the premise of Siegfried Thomas Bok's Law
of Volume Constancy, as disclosed in his 1929 publication on cortical laminae.

"""
from odbAccess import *
from numpy import sum
import math
import csv

#////////////////////// SET-UP
class odb:
   def __init__(self, odbName, part, stepSelection, elementType, elementLayer, setName, 
      upperNodeSetName, lowerNodeSetName, stateVariable, xCoordLower, xCoordUpper, 
      adjustmentLower, adjustmentUpper, adjustmentUpperCoords, connectivityLower, connectivityUpper, frame):
      self.odbName = odbName
      self.odb = openOdb(odbName, readOnly=False)
      self.part = self.odb.rootAssembly[part]
      self.stepSelection = self.odb.steps[stepSelection]
      self.elementType = elementType
      self.elementLayer = self.part.elementSets[elementLayer]
      self.setName = setName
      self.upperNodeSetName = upperNodeSetName
      self.lowerNodeSetName = lowerNodeSetName
      self.stateVariable = stateVariable
      self.xCoordLower = xCoordLower
      self.xCoordUpper = xCoordUpper
      self.adjustmentLower = adjustmentLower
      self.adjustmentUpper = adjustmentUpper 
      self.adjustmentUpperCoords = adjustmentUpperCoords
      self.connectivityLower = connectivityLower
      self.connectivityUpper = connectivityUpper
      self.stateVariableSum = 0
      self.meanLaminaThickness = 0
      self.selectedElementLabels = []
      self.upperNodeList = []
      self.lowerNodeList = []
      self.frame = frame

#///// Function for obtaining element & element-node-connectivity information for determing element bounds of sulci and gyri
   def get_element_info(self):
      upperNodeLeft = []
      lowerNodeLeft = []
      upperNodeRight = []
      lowerNodeRight = [] 
      elementLabels = []
      for element in self.elementLayer.elements:
         elementLabels.append(element.label)
         upperNodeLeft.append(element.connectivity[6])
         lowerNodeLeft.append(element.connectivity[5])
         upperNodeRight.append(element.connectivity[2])
         lowerNodeRight.append(element.connectivity[1])
      outFile = open(self.odbName + '_information_for_gyri_and_sulci_selection.csv', 'w+')
      header = ['Element Label', 'Upper Node Left Number', 'Lower Node Left Number', 'Upper Node Right Number', 'Lower Node Right Number']
      output = zip(elementLabels , upperNodeLeft, lowerNodeLeft, upperNodeRight, lowerNodeRight)
      writer = csv.writer(outFile)
      writer.writerow(header)
      writer.writerows(output)
      outFile.close()

#///// Function for creating element & node sections for sulci and gyri
   def elementSelection(self):
      for element in self.elementLayer.elements:
         if element.label >= self.xCoordLower and element.label <= self.xCoordUpper:
            self.selectedElementLabels.append(element.label)

   def nodeSelection(self):
# lower
      elementBoundLower = self.xCoordLower
      for element in self.elementLayer.elements:
         if element.label >= elementBoundLower and element.label <= (elementBoundLower+1):
            lowerNode = element.connectivity[self.connectivityLower]
            self.lowerNodeList.append(lowerNode)
            print elementBoundLower
            print element.label
            if elementBoundLower < self.xCoordUpper:
               oldbound = elementBoundLower
               elementBoundLower = round(oldbound  + 6,3)
            else:
               break
      print self.lowerNodeList
# Upper
      elementBoundUpper = (self.xCoordLower+5)
      for element in self.elementLayer.elements:
         if element.label >= elementBoundUpper and element.label <= (elementBoundUpper+1):
            upperNode = element.connectivity[self.connectivityUpper]
            self.upperNodeList.append(upperNode)
            print elementBoundUpper
            print element.label
            if elementBoundUpper < self.xCoordUpper:
               oldbound = elementBoundUpper
               elementBoundUpper = round(oldbound  + 6,3)
            else:
               break
      print self.upperNodeList

#///// Function for extracting state variable information from created element set
   def solveStateVariable(self):
      if self.frame == self.stepSelection.frames[0]:
         elementSet = self.part.ElementSetFromElementLabels(name=self.setName, elementLabels=self.selectedElementLabels)
      else:
         elementSet = self.part.elementSets[self.setName]
      frameSelection = self.frame
      stateVariableField = frameSelection.fieldOutputs[self.stateVariable].getSubset(region=elementSet, elementType =self.elementType)
      for j in stateVariableField.values:
         self.stateVariableSum += j.data 

#///// Function for extracting thickness information from created node sets
   def solveThickness(self):
      if self.frame == self.stepSelection.frames[0]:
         upperNodeSet = self.part.NodeSetFromNodeLabels(name=self.upperNodeSetName, nodeLabels=self.upperNodeList)
         lowerNodeSet = self.part.NodeSetFromNodeLabels(name=self.lowerNodeSetName, nodeLabels=self.lowerNodeList)
      else:
         upperNodeSet = self.part.nodeSets[self.upperNodeSetName]
         lowerNodeSet = self.part.nodeSets[self.lowerNodeSetName]
      frameSelection = self.frame
      nodeCoords = frameSelection.fieldOutputs['COORD']
      laminaThickness = []
      for j, k in zip(sorted(nodeCoords.getSubset(region=upperNodeSet).values, reverse=self.adjustmentUpperCoords), nodeCoords.getSubset(region=lowerNodeSet).values):
         print 'x coords'
         print j.data[0], ' - ', k.data[0]
         print 'y coords'
         print j.data[1], ' - ', k.data[1]
         thickness = round(math.sqrt(((j.data[0]-k.data[0]) ** 2)+((j.data[1]-k.data[1]) ** 2)), 8)
         laminaThickness.append(thickness)
      print laminaThickness
      self.meanLaminaThickness = sum(laminaThickness)/len(laminaThickness)

#///// Function for writing state variable and thickness information to file for plotting
def writeOutputToFile(fileName, inputList1):
   if inputList1 == listOfStateVariableList:
      outFile = open(fileName + '_State_Variable_Output.csv', 'w+')
      header = ['Frame', 'Cortical Section', 'Section Type', 'Layer', 'Strain Energy']
      output = zip(frameList, listOfSectionList, listOfSectionTypeList, listOfLayerList, inputList1)
   elif inputList1 == listOfLaminaThicknessList: 
      outFile = open(fileName + '_Lamina_Thickness_Output.csv', 'w+')
      header = ['Frame', 'Cortical Section', 'Section Type', 'Layer', 'Lamina Thickness']
      output = zip(frameList, listOfSectionList, listOfSectionTypeList, listOfLayerList, inputList1)
   writer = csv.writer(outFile)
   writer.writerow(header)
   writer.writerows(output)
   outFile.close()

#///// Classes which initalizes properties of "Example" simulations
class example_simulation:
   def __init__(self):
      self.odbName = 'Rec_6.odb'
      self.lamina = ['CORTEX']
      self.strainEnergy = 'ELSE'
      self.part = 'REC_6'
      self.step = 'pressure'
      self.elementTypeNumber = 'C3D8'
      self.adjustmentLower = [0]
      self.adjustmentUpper = [0]
      self.adjustmentUpperCoords = [True]
      self.connectivityLower = [5]
      self.connectivityUpper = [6]
      self.sectionList = ['SULCI1.1', 'GYRI1.1']
      self.xCoordLowerList = [234.5, 540.5]
      self.xCoordUpperList = [264.5, 570.5]
# This is used to get the state variable information and thickness across all frames
      self.odb = openOdb(self.odbName, readOnly=False)
      self.frames = self.odb.steps[self.step].frames

#////////////////////// SOLVING
def element_info(odbFile):
   analysis = odb(odbFile.odbName, odbFile.part, odbFile.step, odbFile.elementTypeNumber, 'CORTEX', 'temp', 'temp', 'temp', odbFile.strainEnergy, 1, 2, 0, 0, [True], 6, 5, odbFile.frames[-1])
   analysis.get_element_info()

def solve(odbFile):
   global listOfStateVariableList
   global listOfLaminaThicknessList
   global listOfSectionList
   global listOfLayerList
   global listOfSectionTypeList
   global frameList
   listOfStateVariableList = []
   listOfLaminaThicknessList = []
   listOfSectionList = []
   listOfSectionTypeList = []
   listOfLayerList = []
   frameList = []
   for section, j, k in zip(odbFile.sectionList, odbFile.xCoordLowerList, odbFile.xCoordUpperList):
      print section
      for frame in range(len(odbFile.frames)):
         for layer, a, b, c, d, e in zip(odbFile.lamina, odbFile.adjustmentLower, odbFile.adjustmentUpper, odbFile.adjustmentUpperCoords, odbFile.connectivityLower, odbFile.connectivityUpper):
            analysis = odb(odbFile.odbName, odbFile.part, odbFile.step, odbFile.elementTypeNumber, layer, layer + '_' + section, layer + '_' + section + '_UPPER', layer + '_' + section + '_LOWER', odbFile.strainEnergy, j, k, a, b, c, d, e, odbFile.frames[frame])
            analysis.elementSelection()
            analysis.nodeSelection()
# Create section and layer list
            if section[0:10] == 'HORIZONTAL':
               listOfSectionTypeList.append(section[0:10])
            elif section[0:5] == 'SULCI':
               listOfSectionTypeList.append(section[0:5])
            elif section[0:4] == 'GYRI':
               listOfSectionTypeList.append(section[0:4])
            elif section == 'FULL':
               listOfSectionTypeList.append(section)
            listOfSectionList.append(section)
            listOfLayerList.append(layer)
# Solve for state variable and thickness information
            analysis.solveStateVariable()
            listOfStateVariableList.append(analysis.stateVariableSum)
            analysis.solveThickness()
            listOfLaminaThicknessList.append(analysis.meanLaminaThickness)
         frameList.append(frame)
# Write information to file
   writeOutputToFile(odbFile.odbName, listOfStateVariableList)
   writeOutputToFile(odbFile.odbName, listOfLaminaThicknessList)

# Solve for simulations
element_info(example_simulation())
solve(example_simulation())
