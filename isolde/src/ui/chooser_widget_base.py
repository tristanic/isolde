# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'chooser_widget_base.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from Qt import QtCore, QtGui, QtWidgets

class ChooserWindowBase:
    def __init__(self, session, title, chooser_class):
        from chimerax.ui import MainToolWindow
        self.tool_window = tw = MainToolWindow(self)
        parent = tw.ui_area
        from Qt.QtWidgets import QVBoxLayout, QHBoxLayout, QDialogButtonBox,  


class Ui_ChooserWidget(object):
    def setupUi(self, ChooserWidget):
        ChooserWidget.setObjectName("ChooserWidget")
        ChooserWidget.resize(400, 300)
        self.verticalLayout = QtWidgets.QVBoxLayout(ChooserWidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.main_widget = QtWidgets.QWidget(ChooserWidget)
        self.main_widget.setObjectName("main_widget")
        self.verticalLayout.addWidget(self.main_widget)
        self.horizontalFrame = QtWidgets.QFrame(ChooserWidget)
        self.horizontalFrame.setMaximumSize(QtCore.QSize(16777215, 50))
        self.horizontalFrame.setObjectName("horizontalFrame")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.horizontalFrame)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.ok_button = QtWidgets.QPushButton(self.horizontalFrame)
        self.ok_button.setObjectName("ok_button")
        self.horizontalLayout.addWidget(self.ok_button)
        self.cancel_button = QtWidgets.QPushButton(self.horizontalFrame)
        self.cancel_button.setObjectName("cancel_button")
        self.horizontalLayout.addWidget(self.cancel_button)
        self.verticalLayout.addWidget(self.horizontalFrame)

        self.retranslateUi(ChooserWidget)
        QtCore.QMetaObject.connectSlotsByName(ChooserWidget)

    def retranslateUi(self, ChooserWidget):
        _translate = QtCore.QCoreApplication.translate
        ChooserWidget.setWindowTitle(_translate("ChooserWidget", "GroupBox"))
        ChooserWidget.setTitle(_translate("ChooserWidget", "GroupBox"))
        self.ok_button.setText(_translate("ChooserWidget", "OK"))
        self.cancel_button.setText(_translate("ChooserWidget", "Cancel"))

