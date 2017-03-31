# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'molprobity_master.ui'
#
# Created by: PyQt5 UI code generator 5.8
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Frame(object):
    def setupUi(self, Frame):
        Frame.setObjectName("Frame")
        Frame.resize(536, 830)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Frame.sizePolicy().hasHeightForWidth())
        Frame.setSizePolicy(sizePolicy)
        Frame.setMinimumSize(QtCore.QSize(480, 600))
        Frame.setFrameShape(QtWidgets.QFrame.NoFrame)
        Frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(Frame)
        self.verticalLayout_2.setContentsMargins(2, 2, 2, 2)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalFrame = QtWidgets.QFrame(Frame)
        self.horizontalFrame.setMaximumSize(QtCore.QSize(16777215, 50))
        self.horizontalFrame.setObjectName("horizontalFrame")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.horizontalFrame)
        self.horizontalLayout.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(self.horizontalFrame)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.selected_model_combo_box = QtWidgets.QComboBox(self.horizontalFrame)
        self.selected_model_combo_box.setObjectName("selected_model_combo_box")
        self.horizontalLayout.addWidget(self.selected_model_combo_box)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.verticalLayout_2.addWidget(self.horizontalFrame)
        self.main_widget_container = QtWidgets.QVBoxLayout()
        self.main_widget_container.setObjectName("main_widget_container")
        self.verticalLayout_2.addLayout(self.main_widget_container)

        self.retranslateUi(Frame)
        QtCore.QMetaObject.connectSlotsByName(Frame)

    def retranslateUi(self, Frame):
        _translate = QtCore.QCoreApplication.translate
        Frame.setWindowTitle(_translate("Frame", "Frame"))
        self.label.setText(_translate("Frame", "Selected model: "))

