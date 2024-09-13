
#QtCore is the nonGUI stuff.


# Import required modules
import sys, time, os

from PySide6.QtWidgets import QMenu, QApplication, QStatusBar, QMessageBox, QFileDialog
from PySide6.QtCore import QSettings
from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.QtCore import QUrl
from PySide6.QtGui import QDesktopServices, QIcon, QAction, QKeySequence
from PySide6.QtWidgets import QLabel


#Import model
from rdeditor.molEditWidget import MolEditWidget
from rdeditor.molViewWidget import MolWidget
from rdeditor.ptable_widget import PTable

from rdkit import Chem


from sddata import SDdata



#The Main browser windows
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self,  fileName=None):
        super(MainWindow,self).__init__()
        self.fileName = fileName
        self.filter = "SD files (*.sdf *.sd)"
        self.model = SDdata()
        #Setup the user interface
        #self.initUI()

        #If we get a filename, load it into the model
        if self.fileName != None:
            self.model.loadSDfile(fileName)


        #Set Window properties
        self.setWindowTitle("A Simple SD file browser")
        self.setWindowIcon(QIcon('Peptide.png'))
        self.setGeometry(100, 100, 200, 150)
        #Set Central Widget
        # self.center = QtSvg.QSvgWidget()
        # self.center.setFixedSize(350,350)
        self.molviewer = MolWidget()
        self.molviewer.setFixedSize(600, 600)
        self.setCentralWidget(self.molviewer)
        #Setup the statusbar
        self.myStatusBar = QStatusBar()
        #A permanent widget is right aligned
        self.molcounter = QLabel("-/-")
        self.myStatusBar.addPermanentWidget(self.molcounter, 0)
        self.setStatusBar(self.myStatusBar)
        self.myStatusBar.showMessage('Ready', 10000)
        self.define_actions()
        self.define_menus()

        #Connect model signals to UI slots
        #Update central widget if the selected molecule changes
        self.model.selectedChanged.connect(self.update_mol)
        #Update the permanent widget in the status bar, if status changes
        self.model.statusChanged.connect(self.molcounter.setText)
        #Finally! Show the UI!
        self.update_mol()
        self.show()


    def update_mol(self):
        self.molviewer.mol = self.model.get_selected_mol()


    #Open a new file
    def openFile(self):
        self.fileName, self.filter = QFileDialog.getOpenFileName(self, filter=self.filter)
        self.model.loadSDfile(str(self.fileName))
    #Increment the selected mol with 1
    def nextMol(self):
        #Increment the selected molecule in the model by 1
        self.model.setSelected(self.model.selected + 1)
    #Decrement the selected mol with 1
    def prevMol(self):
        #Decrement the selected molecule in the model by 1
        self.model.setSelected(self.model.selected - 1)


    # Menus
    def define_menus(self):
        #Setup the menu
        self.fileMenu = self.menuBar().addMenu("&File")
        self.helpMenu = self.menuBar().addMenu("&Help")
        #Setup the Toolbar
        self.mainToolBar = self.addToolBar('Main')
        #Populate the Menu with Actions
        self.fileMenu.addAction(self.openAction)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAction)
        self.helpMenu.addAction(self.aboutAction)
        self.helpMenu.addSeparator()
        self.helpMenu.addAction(self.aboutQtAction)
        #Populate the Toolbar with actions.
        self.mainToolBar.addAction(self.openAction)
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.prevAction)
        self.mainToolBar.addAction(self.nextAction)
        # self.mainToolBar.addSeparator()
        # self.mainToolBar.addAction(self.molblockAction)


    # ACTIONS
    def define_actions(self):
        self.openAction = QAction( QIcon('Open Folder.png'), 'O&pen',
                            self, shortcut=QKeySequence.Open,
                            statusTip="Open an SD file",
                            triggered=self.openFile)
        
        # self.molblockAction = QAction( QIcon('Page Overview 3.png'), 'V&iew MolBlock',
        #                           self, shortcut="Ctrl+M",
        #                           statusTip="View MolBlock",
        #                           triggered=self.viewMolBlock)
        self.exitAction = QAction( QIcon('Exit.png'), 'E&xit',
                                   self, shortcut="Ctrl+Q",
                                   statusTip="Close the Application",
                                   triggered=self.exitFile)
        self.prevAction = QAction( QIcon('Left.png'),'Previous', self,
                                   shortcut=QKeySequence.MoveToPreviousChar,
                                   statusTip="Previous molecule",
                                   triggered=self.prevMol)
        self.nextAction = QAction( QIcon('Right.png'),'Next', self,
                                   shortcut=QKeySequence.MoveToNextChar,
                                   statusTip="Next molecule",
                                   triggered=self.nextMol)
        self.aboutAction = QAction( QIcon('Info.png'), 'A&;bout',
                                    self, statusTip="Displays info about SDbrowser",
                                   triggered=self.aboutHelp)
        self.aboutQtAction = QAction("About &Qt", self,
                                statusTip="Qt library About box",
                                triggered=qApp.aboutQt)
        

    def exitFile(self):
        response = self.msgApp("Confirmation", "This will quit the application. Do you want to Continue?")
        if response == "Y":
            self.ptable.close()
            exit(0)
        else:
            self.editor.logger.debug("Abort closing")

    def aboutHelp(self):
        QMessageBox.about(self, "A Basic SD browser",
                "Hackathon work on a GUI for sdfile unsanitazion fixer")
        

    # Function to show Diaglog box with provided Title and Message
    def msgApp(self, title, msg):
        userInfo = QMessageBox.question(self, title, msg)
        if userInfo == QMessageBox.Yes:
            return "Y"
        if userInfo == QMessageBox.No:
            return "N"
        self.close()


if __name__ == '__main__':
    # Exception Handling
    try:
        sdBrowser = QApplication(sys.argv)
        #Load with file if provided
        if len(sys.argv) > 1:
            mainWindow = MainWindow(fileName = sys.argv[1])
        else:
            mainWindow = MainWindow()
        sdBrowser.exec()
        sys.exit(0)
    #Basic Exception handling
    except NameError:
        print("Name Error:", sys.exc_info()[1])
    except SystemExit:
        print("Closing")
    except Exception:
        print(sys.exc_info()[1])