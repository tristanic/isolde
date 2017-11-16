# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'src/IsoldeWindow.ui'
#
# Created by: PyQt5 UI code generator 5.9.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Isolde_MainWin(object):
    def setupUi(self, Isolde_MainWin):
        Isolde_MainWin.setObjectName("Isolde_MainWin")
        Isolde_MainWin.resize(574, 863)
        self.centralwidget = QtWidgets.QWidget(Isolde_MainWin)
        self.centralwidget.setObjectName("centralwidget")
        Isolde_MainWin.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(Isolde_MainWin)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 574, 28))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        Isolde_MainWin.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(Isolde_MainWin)
        self.statusbar.setObjectName("statusbar")
        Isolde_MainWin.setStatusBar(self.statusbar)
        self.actionSave_current_model_as = QtWidgets.QAction(Isolde_MainWin)
        self.actionSave_current_model_as.setObjectName("actionSave_current_model_as")
        self.actionLoad_example_data = QtWidgets.QAction(Isolde_MainWin)
        self.actionLoad_example_data.setObjectName("actionLoad_example_data")
        self.actionLoad_after_example_data = QtWidgets.QAction(Isolde_MainWin)
        self.actionLoad_after_example_data.setObjectName("actionLoad_after_example_data")
        self.menuFile.addAction(self.actionSave_current_model_as)
        self.menuFile.addAction(self.actionLoad_example_data)
        self.menuFile.addAction(self.actionLoad_after_example_data)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(Isolde_MainWin)
        QtCore.QMetaObject.connectSlotsByName(Isolde_MainWin)

    def retranslateUi(self, Isolde_MainWin):
        _translate = QtCore.QCoreApplication.translate
        Isolde_MainWin.setWindowTitle(_translate("Isolde_MainWin", "MainWindow"))
        self.menuFile.setTitle(_translate("Isolde_MainWin", "File"))
        self.actionSave_current_model_as.setText(_translate("Isolde_MainWin", "Save current model as..."))
        self.actionLoad_example_data.setText(_translate("Isolde_MainWin", "Load \"before\" example data"))
        self.actionLoad_after_example_data.setText(_translate("Isolde_MainWin", "Load \"after\" example data"))

from . import resources_rc
