# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ExcelWindow2.ui'
#
# Created by: PyQt5 UI code generator 5.15.6
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(711, 293)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(60, 40, 66, 16))
        self.label.setObjectName("label")
        self.excelPath = QtWidgets.QLineEdit(self.centralwidget)
        self.excelPath.setGeometry(QtCore.QRect(140, 40, 391, 20))
        self.excelPath.setObjectName("excelPath")
        self.chooseExcel = QtWidgets.QPushButton(self.centralwidget)
        self.chooseExcel.setGeometry(QtCore.QRect(560, 40, 61, 23))
        self.chooseExcel.setObjectName("chooseExcel")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(60, 90, 60, 16))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(130, 90, 48, 16))
        self.label_3.setObjectName("label_3")
        self.kLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.kLineEdit.setGeometry(QtCore.QRect(220, 90, 61, 20))
        self.kLineEdit.setObjectName("kLineEdit")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(300, 90, 84, 20))
        self.label_4.setObjectName("label_4")
        self.stepLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.stepLineEdit.setGeometry(QtCore.QRect(390, 90, 61, 20))
        self.stepLineEdit.setObjectName("stepLineEdit")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(470, 90, 84, 16))
        self.label_5.setObjectName("label_5")
        self.zMaxLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.zMaxLineEdit.setGeometry(QtCore.QRect(560, 90, 61, 20))
        self.zMaxLineEdit.setObjectName("zMaxLineEdit")
        self.start = QtWidgets.QPushButton(self.centralwidget)
        self.start.setGeometry(QtCore.QRect(560, 170, 61, 23))
        self.start.setObjectName("start")
        self.xBackLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.xBackLineEdit.setGeometry(QtCore.QRect(220, 130, 61, 20))
        self.xBackLineEdit.setObjectName("xBackLineEdit")
        self.yBackLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.yBackLineEdit.setGeometry(QtCore.QRect(390, 130, 61, 20))
        self.yBackLineEdit.setObjectName("yBackLineEdit")
        self.zBackLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.zBackLineEdit.setGeometry(QtCore.QRect(560, 130, 61, 20))
        self.zBackLineEdit.setObjectName("zBackLineEdit")
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(130, 130, 91, 20))
        self.label_6.setObjectName("label_6")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(300, 130, 91, 20))
        self.label_7.setObjectName("label_7")
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setGeometry(QtCore.QRect(470, 130, 91, 20))
        self.label_8.setObjectName("label_8")
        self.label_9 = QtWidgets.QLabel(self.centralwidget)
        self.label_9.setGeometry(QtCore.QRect(130, 170, 91, 20))
        self.label_9.setObjectName("label_9")
        self.dLineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.dLineEdit.setGeometry(QtCore.QRect(220, 170, 61, 20))
        self.dLineEdit.setObjectName("dLineEdit")
        self.label_10 = QtWidgets.QLabel(self.centralwidget)
        self.label_10.setGeometry(QtCore.QRect(300, 170, 71, 20))
        self.label_10.setObjectName("label_10")
        self.d2LineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.d2LineEdit.setGeometry(QtCore.QRect(390, 170, 61, 20))
        self.d2LineEdit.setObjectName("d2LineEdit")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 711, 23))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.chooseExcel.clicked.connect(
            MainWindow.chooseExcel)  # type: ignore
        self.start.clicked.connect(MainWindow.start)  # type: ignore
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "gcode生成助手"))
        self.label.setText(_translate("MainWindow", "dxf路径："))
        self.chooseExcel.setText(_translate("MainWindow", "选择"))
        self.label_2.setText(_translate("MainWindow", "参数设置："))
        self.label_3.setText(_translate("MainWindow", "倍率：k="))
        self.label_4.setText(_translate("MainWindow", "z增加量：step="))
        self.label_5.setText(_translate("MainWindow", "z最大值：zMax="))
        self.start.setText(_translate("MainWindow", "生成"))
        self.label_6.setText(_translate("MainWindow", "x退刀量：xBack="))
        self.label_7.setText(_translate("MainWindow", "y退刀量：yBack="))
        self.label_8.setText(_translate("MainWindow", "z退刀量：zBack="))
        self.label_9.setText(_translate("MainWindow", "重边偏移量：d="))
        self.label_10.setText(_translate("MainWindow", "延长量：d2="))