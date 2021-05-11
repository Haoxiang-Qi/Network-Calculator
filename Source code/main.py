import sys
from PyQt5.QtWidgets import QApplication, QMainWindow
from tool_new import *
import  multiprocessing
from PyQt5 import QtGui
from PyQt5.QtGui import *

class MyWindow(QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super(MyWindow, self).__init__(parent)
        self.setupUi(self)
    def open(self):
        self.show()
class New_win1(QMainWindow,Ui_MainWindow1):
    def __init__(self, parent=None):
        super(New_win1, self).__init__(parent)
        self.setupUi(self)
    def open(self):
        self.show()
class New_win2(QMainWindow,Ui_MainWindow2):
    def __init__(self, parent=None):
        super(New_win2, self).__init__(parent)
        self.setupUi(self)
    def open(self):
        self.show()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    myWin = MyWindow()
    myWin.setWindowTitle("Network calculator")
    win1 = New_win1()
    win2 = New_win2()
    myWin.actionMaximum_connection_value_analysis.triggered.connect(win1.open)
    myWin.action_window3.triggered.connect(win2.open)
    win1.actionwindow3.triggered.connect(win2.open)
    win1.actionGene_network_connect.triggered.connect(myWin.open)
    win2.actionMaxium.triggered.connect(win1.open)
    win2.actionGene_network_connect.triggered.connect(myWin.open)
    myWin.show()
    sys.exit(app.exec_())