#QtCore is the nonGUI stuff.
from PySide6.QtCore import QObject, Slot
from PySide6.QtCore import Signal
from PySide6 import QtWidgets
#RDKit stuff
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import PandasTools

import sys

import logging
logger = logging.getLogger()
#The model holding the SDfile data
class SDdata(QObject): #Inherit from QObject so that Signals can be emitted
    def __init__(self):
        super().__init__() #Init the super QObject class so that it works with QT stuff etc.
        self._selected = 0
        self._status = "Ready"
        #self.length = 0

    selectedChanged = Signal(int, name = 'selectedChanged')
    statusChanged = Signal(int, name = "statusChanged")

    @property
    def selected(self):
        return self._selected
    #This enables us to do model.selected = 2
    @selected.setter
    def selected(self, value):
        self.setSelected(value)
    #this is more easy to use from with PyQt signals
    def setSelected(self, selected):
        #Prevent setting a selected that doesn't exist
        if selected < 0: selected = 0
        if selected > len(self) -1: selected = len(self) -1
        #Only set the selected if its changed, we get round tripping otherwise
        if selected != self._selected:
            self._selected = selected
            logger.debug("in model: selected set for ", selected)
            #Emit the signal that selected has changed
            self.selectedChanged.emit(self._selected)

    #Decorate the function setCounter
    @Slot()
    def setCounter(self):
        self.setStatus('%s/%s'%(self._selected + 1, len(self)))

        #Set the counter when the selection is changed (non GUI signal slot example)
        self.selectedChanged.connect(self.setCounter)

    def __len__(self):
        return len(self.sddata)

    #TODO: maybe rather load into a dataframe??
    def loadSDfile(self, filename):
        self.filename = filename
        #self.SDMolSupplier = Chem.SDMolSupplier(filename, sanitize=False, strictParsing=False)

        self.sddata = PandasTools.LoadSDF(filename, removeHs=False, strictParsing=False, sanitize=False)
        
        if self.selected == 0:
            self.selectedChanged.emit(self._selected)
        else:
            self.setSelected(0)


    def get_selected_mol(self):
        return self.sddata.ROMol.iloc[self._selected]