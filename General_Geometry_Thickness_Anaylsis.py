"""
Written by: Jack Consolini
Last updated: 12/06/2022


This python script is used to extract laminar volume and thickness information from 
2D multi-layered simulations of brain (cortex) folding. This script was used for the
conference poster presentation: "Investigation of law of volume constancy in cortical 
lamina", and offered insight into my recently published review paper, "Bok's equi-volume
principle: translation, historical context, and a modern perspective".

In order to extract laminar volume and thickness information from 2D layered output database 
(odb) files, the functions accomplish the following:

"""

#////////////////////// IMPORTS
# Allow access to odb files from 2D layered simulations
from odbAccess import *
from numpy import sum
import math
import csv

#////////////////////// SET-UP
# Create class to initialize all finite element node and element information saved 
# within the odb file.
class odb:
   def __init__(self, odbName, part, stepSelection, elementType, elementLayer, setName, 
      upperNodeSetName, lowerNodeSetName, stateVariable, xCoordLower, xCoordUpper, 
      adjustmentLower, adjustmentUpper, adjustmentUpperCoords, connectivityLower, 
      connectivityUpper, frame):
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

#////////////////////// NODES ASSOCIATED TO FINITE ELEMENTS
# Function for obtaining the nodes that are connected to the finite elements. 
# The volume and thickness of these elements will be solved for each layer.
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
# Write element-node connectivity to .csv for quality checking.
      outFile = open(self.odbName + '_information_for_gyri_and_sulci_selection.csv', 'w+')
      header = ['Element Label', 'Upper Node Left Number', 'Lower Node Left Number', 
      'Upper Node Right Number', 'Lower Node Right Number']
      output = zip(elementLabels , upperNodeLeft, lowerNodeLeft, upperNodeRight, lowerNodeRight)
      writer = csv.writer(outFile)
      writer.writerow(header)
      writer.writerows(output)
      outFile.close()

#////////////////////// SELECT SPECIFIC ELEMENTS WITHIN A LAYER
# Function for selecting specific elements within a layer, to analyze the volume and 
# thickness at the top and bottom of folded sections (straight sections were excluded 
# from the analysis).
   def elementSelection(self):
      for element in self.elementLayer.elements:
         if element.label >= self.xCoordLower and element.label <= self.xCoordUpper:
            self.selectedElementLabels.append(element.label)

#////////////////////// CREATE NODE SETS FOR THICKNESS CALCULATION
# Function for creating sets of nodes along the bottom and top of the chosen elements, of which
# the distance between these nodes will be used to determine the thickness of the elements.
   def nodeSelection(self):
# Lower node set
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
# Upper node set
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

#////////////////////// DETERMINE THICKNESS OF SELECTED ELEMENTS
# Function for extracting element thickness using the created node sets.
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
      for j, k in zip(sorted(nodeCoords.getSubset(region=upperNodeSet).values, reverse=self.adjustmentUpperCoords), 
         nodeCoords.getSubset(region=lowerNodeSet).values):
         print 'x coords'
         print j.data[0], ' - ', k.data[0]
         print 'y coords'
         print j.data[1], ' - ', k.data[1]
         thickness = round(math.sqrt(((j.data[0]-k.data[0]) ** 2)+((j.data[1]-k.data[1]) ** 2)), 8)
         laminaThickness.append(thickness)
      print laminaThickness
      self.meanLaminaThickness = sum(laminaThickness)/len(laminaThickness)

#////////////////////// EXTRACT VOLUME INFORMATION FROM SELECTED ELEMENTS
# Function for extracting volume (which is a state variable) from the set of 
# selected elements.
   def solveStateVariable(self):
      if self.frame == self.stepSelection.frames[0]:
         elementSet = self.part.ElementSetFromElementLabels(name=self.setName, 
            elementLabels=self.selectedElementLabels)
      else:
         elementSet = self.part.elementSets[self.setName]
      frameSelection = self.frame
      stateVariableField = frameSelection.fieldOutputs[self.stateVariable].getSubset(region=elementSet, 
         elementType =self.elementType)
      for j in stateVariableField.values:
         self.stateVariableSum += j.data 

#////////////////////// WRITE VOLUME AND THICKNESS DATA TO FILE
# Function for writing volume and thickness information to file for plotting.
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

#////////////////////// EXAMPLE ABAQUS SIMULATION
# This class initalizes an odb file to serve as an example of how the code works. The .odb 
# and an Abaqus liscense is required for the code to succeed. This class shows what user 
# inputs are required for the code to run; such as the xCoordLowerList and xCoordUpperList, 
# which tell the script which elements to extract data from.
class example_simulation:
   def __init__(self):
      self.odbName = '2D_6L_SR1_A0_Je_v2.odb'
      self.lamina = ['BOX-1_GRAY1', 'BOX-1_GRAY2', 'BOX-1_GRAY3', 'BOX-1_GRAY4', 
      'BOX-1_GRAY5', 'BOX-1_GRAY6']
      self.adjustmentLower = [-2, -1, -1, -1, -1, 0]
      self.adjustmentUpper = [0, 1, -1, -1, -1, 0]
      self.adjustmentUpperCoords = [True, True, False, False, False, False]
      self.connectivityLower = [0, 3, 3, 3, 3, 0]
      self.connectivityUpper = [3, 0, 0, 0, 0, 3]
      self.sectionList = ['HORIZONTAL1', 'SULCI1', 'HORIZONTAL2', 'GYRI1', 
      'HORIZONTAL3', 'SULCI2', 'HORIZONTAL4', 'GYRI2', 'HORIZONTAL5', 'SULCI3', 
      'HORIZONTAL6', 'GYRI3', 'HORIZONTAL7', 'SULCI4', 'HORIZONTAL8', 'GYRI4', 
      'HORIZONTAL9', 'SULCI5', 'HORIZONTAL10', 'GYRI5', 'HORIZONTAL11', 'SULCI6', 
      'HORIZONTAL12', 'GYRI6', 'HORIZONTAL13', 'SULCI7', 'HORIZONTAL14' ]
      self.xCoordLowerList = [-0.0272, -0.0260, -0.0248, -0.0216, -0.0184, -0.0172, -0.0160, -0.0132, -0.0099, -0.0088, -0.0076, -0.0048, -0.0019, -0.0009, 0.0004, 0.0032, 0.0060, 0.0072, 0.0084, 0.0112, 0.0144, 0.0156, 0.0168, 0.0200, 0.0236, 0.0248, 0.0260]
      self.xCoordUpperList = [-0.0260, -0.0248, -0.0236, -0.0204, -0.0172, -0.0160, -0.0148, -0.0120, -0.0088, -0.0076, -0.0064, -0.0036, -0.0009, 0.0004, 0.0016, 0.0044, 0.0072, 0.0084, 0.0096, 0.0124, 0.0156, 0.0168, 0.0180, 0.0212, 0.0248, 0.0260, 0.0272]
# This is used to get the state variable information and thickness across all frames
      self.odb = openOdb(self.odbName, readOnly=False)
      self.frames = self.odb.steps[self.step].frames

#////////////////////// SOLVING
# This function first obtains the information of what nodes are associated to what elements.
def element_info(odbFile):
   analysis = odb(odbFile.odbName, odbFile.part, odbFile.step, odbFile.elementTypeNumber, 'CORTEX', 'temp', 
      'temp', 'temp', odbFile.strainEnergy, 1, 2, 0, 0, [True], 6, 5, odbFile.frames[-1])
   analysis.get_element_info()
# This function solves for the volume and thickness in the specified sections, which were 
# normally made of up of 4-5 finite elements, in each layer. Once the information for each 
# section in each layer is determined, the relative volume (section volume per layer over 
# section volume for all layers) and the relative thickness (section thickness per layer 
# over section thickness for all layers) is calculated. Finally, this information is then 
# saved to a .csv file for plotting with an additional scipt.
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
         for layer, a, b, c, d, e in zip(odbFile.lamina, odbFile.adjustmentLower, odbFile.adjustmentUpper, 
            odbFile.adjustmentUpperCoords, odbFile.connectivityLower, odbFile.connectivityUpper):
            analysis = odb(odbFile.odbName, odbFile.part, odbFile.step, odbFile.elementTypeNumber, layer, 
               layer + '_' + section, layer + '_' + section + '_UPPER', layer + '_' + section + '_LOWER', 
               odbFile.strainEnergy, j, k, a, b, c, d, e, odbFile.frames[frame])
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
# Solve for volume and thickness information
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
