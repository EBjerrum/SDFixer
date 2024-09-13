#QtCore is the nonGUI stuff.
from PySide6.QtCore import QObject, Slot
from PySide6.QtCore import Signal
from PySide6 import QtWidgets
#RDKit stuff
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
#The model holding the SDfile data
class SDdata(QObject): #Inherit from QObject so that Signals can be emitted
    def __init__(self):
        super().__init__() #Init the super QObject class so that it works with QT stuff etc.
        self._selected = 0
        self._status = "Ready"
        self.length = 0

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
        if selected > self.length -1: selected = self.length -1
        #Only set the selected if its changed, we get round tripping otherwise
        if selected != self._selected:
            self._selected = selected
            print("in model: selected set for ", selected)
            #Emit the signal that selected has changed
            self.selectedChanged.emit(self._selected)

    #Decorate the function setCounter
    @Slot()
    def setCounter(self):
        self.setStatus('%s/%s'%(self._selected + 1, self.length))

        #Set the counter when the selection is changed (non GUI signal slot example)
        self.selectedChanged.connect(self.setCounter)


    #TODO: maybe rather load into a dataframe??
    def loadSDfile(self, filename):
        self.filename = filename
        self.SDMolSupplier = Chem.SDMolSupplier(filename)
        self.length = len(self.SDMolSupplier)
        if self.selected == 0:
            self.selectedChanged.emit(self._selected)
        else:
            self.setSelected(0)


    def get_selected_mol(self):
        return self.SDMolSupplier[self._selected]

    # #Better rendering with SVG
    # def getMolSvg(self, kekulize=True, calc2Dcoords=True):
    #     mol = self.SDMolSupplier[self._selected]
    #     mc = Chem.Mol(mol.ToBinary())
    #     if kekulize:
    #         try:
    #             Chem.Kekulize(mc)
    #         except:
    #             mc = Chem.Mol(mol.ToBinary())
    #     if not mc.GetNumConformers() or calc2Dcoords:
    #         rdDepictor.Compute2DCoords(mc)
    #     drawer = rdMolDraw2D.MolDraw2DSVG(300,300)
    #     drawer.DrawMolecule(mc)
    #     drawer.FinishDrawing()
    #     svg = drawer.GetDrawingText().replace('svg:','')
    #     return svg
    # def getMolBlock(self):
    #     print("Getting MolBlock") #TODO add logger
    #     return self.SDMolSupplier.GetItemText(self.selected)



if __name__ == '__main__':
    # Exception Handling
    try:
        sdBrowser = QtWidgets.QApplication(sys.argv)
        #Load with file if provided
        if len(sys.argv) > 1:
            mainWindow = MainWindow(fileName = sys.argv[1])
        else:
            mainWindow = MainWindow()
        sdBrowser.exec_()
        sys.exit(0)
    #Basic Exception handling
    except NameError:
        print("Name Error:", sys.exc_info()[1])
    except SystemExit:
        print("Closing")
    except Exception:
        print(sys.exc_info()[1])