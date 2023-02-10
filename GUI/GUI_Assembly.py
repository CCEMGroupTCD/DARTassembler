import sys
from PyQt5.QtWidgets import QMainWindow, QTableWidgetItem
from PyQt5 import QtCore, QtWidgets, QtGui, QtSvg
import json
from src03_Assembly_Cian.main import Assembly

global assembly_input_list
assembly_input_list = []


class Window2(QMainWindow):  # This window corresponds to the Create Batch window
    def setupUi(self, Batch_Widget):  # Here we set up the widget
        self.setFixedSize(800, 600)
        self.border_curvature = "10px"
        self.option_button_height = 30  # 30
        self.option_button_widht = 140
        self.border_width = "3px"
        self.button_spacing = 23
        self.metal_selection_section_gap = 125
        self.viewer_height = 150
        self.viewer_width = 300
        self.left_margin = 40
        self.top_margin = 55
        self.tmp_Dic = {}

        Batch_Widget.setObjectName("Batch_Widget")
        Batch_Widget.resize(800, 600)
        Batch_Widget.setStyleSheet("background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 rgba(5, 105, 185, 255), stop:1 rgba(0, 170, 205, 255));")
        self.Save = QtWidgets.QPushButton(Batch_Widget)
        self.Save.setGeometry(QtCore.QRect(self.left_margin, 518, (2 * self.option_button_widht) + self.button_spacing, self.option_button_height))
        self.Save.setObjectName("Save")
        self.Save.setStyleSheet("QPushButton"
                                "{"
                                "background-color : rgba(255, 210, 0, 255); border-radius: 15px; color: rgba(5, 105, 185, 255); font-size: 20pt"
                                "}"
                                "QPushButton::pressed"
                                "{"
                                "background-color : transparent; color: rgba(255, 210, 0, 255); border-style: outset; border-width: 2px; border-radius: 15px; border-color: rgba(255, 210, 0, 255)"
                                "}")
        self.Save.clicked.connect(self.close_assembly_widget)

        self.comboBox_2 = QtWidgets.QComboBox(Batch_Widget)
        self.comboBox_2.setGeometry(QtCore.QRect(self.left_margin, self.top_margin, round(self.option_button_widht * 0.6), self.option_button_height))
        self.comboBox_2.setObjectName("comboBox_2")
        self.comboBox_2.setStyleSheet("QComboBox::drop-down { "
                                      "subcontrol-origin: padding; "
                                      "subcontrol-position: top right; "
                                      "width: 10px; "
                                      "border-left-width: 5px; "
                                      "border-left-color: rgba(255, 210, 0, 255); "
                                      "border-left-style: solid; "
                                      "border-top-right-radius: 3px; "
                                      "border-bottom-right-radius: 3px; "
                                      "}"
                                      "QComboBox"
                                      "{"
                                      "background-color : white; border-radius: " + self.border_curvature + "; color: rgba(5, 105, 185, 255); "
                                                                                                            "}"
                                                                                                            "QListView"
                                                                                                            "{"
                                                                                                            "background-color: white; selection-color: green"
                                                                                                            "}"
                                                                                                            "QComboBox::item:selected { background-color: rgba(255, 210, 0, 255); color: white }"
                                      )
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_3 = QtWidgets.QComboBox(Batch_Widget)
        self.comboBox_3.setGeometry(QtCore.QRect(round(self.left_margin + (self.option_button_widht * 0.6) + 24), self.top_margin, int(self.option_button_widht * 0.6), self.option_button_height))
        self.comboBox_3.setObjectName("comboBox_3")
        self.comboBox_3.setStyleSheet("QComboBox::drop-down { "
                                      "subcontrol-origin: padding; "
                                      "subcontrol-position: top right; "
                                      "width: 10px; "
                                      "border-left-width: 5px; "
                                      "border-left-color: rgba(255, 210, 0, 255); "
                                      "border-left-style: solid; "
                                      "border-top-right-radius: 3px; "
                                      "border-bottom-right-radius: 3px; "
                                      "}"
                                      "QComboBox"
                                      "{"
                                      "background-color : white; border-radius: " + self.border_curvature + "; color: rgba(5, 105, 185, 255); "
                                                                                                            "}"
                                                                                                            "QListView"
                                                                                                            "{"
                                                                                                            "background-color: white; selection-color: green"
                                                                                                            "}"
                                                                                                            "QComboBox::item:selected { background-color: rgba(255, 210, 0, 255); color: white }"
                                      )
        self.comboBox_3.addItem("")
        self.comboBox_3.addItem("")
        self.comboBox_3.addItem("")
        self.comboBox_3.addItem("")
        self.comboBox_3.addItem("")
        self.comboBox_3.addItem("")
        self.comboBox_3.addItem("")
        self.comboBox_3.addItem("")
        self.comboBox_3.addItem("")
        self.comboBox = QtWidgets.QComboBox(Batch_Widget)
        self.comboBox.setGeometry(QtCore.QRect(round(self.left_margin + (self.option_button_widht * 1.2) + 24 + 24), self.top_margin, int(self.option_button_widht * 0.6), self.option_button_height))
        self.comboBox.setObjectName("comboBox")
        self.comboBox.setStyleSheet("QComboBox::drop-down { "
                                    "subcontrol-origin: padding; "
                                    "subcontrol-position: top right; "
                                    "width: 10px; "
                                    "border-left-width: 5px; "
                                    "border-left-color: rgba(255, 210, 0, 255); "
                                    "border-left-style: solid; "
                                    "border-top-right-radius: 3px; "
                                    "border-bottom-right-radius: 3px; "
                                    "}"
                                    "QComboBox"
                                    "{"
                                    "background-color : white; border-radius: " + self.border_curvature + "; color: rgba(5, 105, 185, 255); "
                                                                                                          "}"
                                                                                                          "QListView"
                                                                                                          "{"
                                                                                                          "background-color: white; selection-color: green"
                                                                                                          "}"
                                                                                                          "QComboBox::item:selected { background-color: rgba(255, 210, 0, 255); color: white }"
                                    )
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.pushButton = QtWidgets.QPushButton(Batch_Widget)
        self.pushButton.setGeometry(QtCore.QRect(self.left_margin, self.top_margin + self.option_button_height + self.button_spacing, self.option_button_widht, self.option_button_height))
        self.pushButton.setObjectName("pushButton")
        self.pushButton.setStyleSheet("QPushButton"
                                      "{"
                                      "background-color : white; border-radius: 10px; color: rgba(5, 105, 185, 255)"
                                      "}"
                                      "QPushButton::pressed"
                                      "{"
                                      "background-color : transparent; color: white; border-style: outset; border-width: 2px; border-radius: " + self.border_curvature + "; border-color: white"
                                                                                                                                                                         "}")
        self.pushButton_2 = QtWidgets.QPushButton(Batch_Widget)
        self.pushButton_2.setGeometry(
            QtCore.QRect(self.left_margin + self.option_button_widht + self.button_spacing, self.top_margin + self.option_button_height + self.button_spacing, self.option_button_widht,
                         self.option_button_height))
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.setStyleSheet("QPushButton"
                                        "{"
                                        "background-color : white; border-radius: 10px; color: rgba(5, 105, 185, 255)"
                                        "}"
                                        "QPushButton::pressed"
                                        "{"
                                        "background-color : transparent; color: white; border-style: outset; border-width: 2px; border-radius: " + self.border_curvature + "; border-color: white"
                                                                                                                                                                           "}")
        self.frame = QtWidgets.QFrame(Batch_Widget)
        self.frame.setGeometry(QtCore.QRect(30, 40, 281, 141))
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.frame.setStyleSheet("background-color : transparent;")
        self.frame_2 = QtWidgets.QFrame(Batch_Widget)
        self.frame_2.setGeometry(QtCore.QRect(30, 210, 281, 141))
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.frame_2.setStyleSheet("background-color : transparent;")
        self.pushButton_3 = QtWidgets.QPushButton(Batch_Widget)
        self.pushButton_3.setGeometry(
            QtCore.QRect(self.left_margin + self.option_button_widht + self.button_spacing, self.top_margin + self.metal_selection_section_gap + self.option_button_height + self.button_spacing,
                         self.option_button_widht, self.option_button_height))
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_3.setShortcut("Del")
        self.pushButton_3.setStyleSheet("QPushButton"
                                        "{"
                                        "background-color : white; border-radius: 10px; color: rgba(5, 105, 185, 255)"
                                        "}"
                                        "QPushButton::pressed"
                                        "{"
                                        "background-color : transparent; color: white; border-style: outset; border-width: 2px; border-radius: " + self.border_curvature + "; border-color: white"
                                                                                                                                                                           "}")
        self.pushButton_4 = QtWidgets.QPushButton(Batch_Widget)
        self.pushButton_4.setGeometry(QtCore.QRect(self.left_margin + self.option_button_widht + self.button_spacing,
                                                   self.top_margin + self.metal_selection_section_gap + (2 * self.option_button_height) + (2 * self.button_spacing), self.option_button_widht,
                                                   self.option_button_height))
        self.pushButton_4.setObjectName("pushButton_4")
        self.pushButton_4.setStyleSheet("QPushButton"
                                        "{"
                                        "background-color : white; border-radius: 10px; color: rgba(5, 105, 185, 255)"
                                        "}"
                                        "QPushButton::pressed"
                                        "{"
                                        "background-color : transparent; color: white; border-style: outset; border-width: 2px; border-radius: " + self.border_curvature + "; border-color: white"
                                                                                                                                                                           "}")
        self.comboBox_4 = QtWidgets.QComboBox(Batch_Widget)
        self.comboBox_4.setGeometry(
            QtCore.QRect(self.left_margin + self.option_button_widht + self.button_spacing, self.top_margin + self.metal_selection_section_gap, self.option_button_widht,
                         self.option_button_height))
        self.comboBox_4.setObjectName("comboBox_4")
        self.comboBox_4.setStyleSheet("QComboBox::drop-down { "
                                      "subcontrol-origin: padding; "
                                      "subcontrol-position: top right; "
                                      "width: 10px; "
                                      "border-left-width: 5px; "
                                      "border-left-color: rgba(255, 210, 0, 255); "
                                      "border-left-style: solid; "
                                      "border-top-right-radius: 3px; "
                                      "border-bottom-right-radius: 3px; "
                                      "}"
                                      "QComboBox"
                                      "{"
                                      "background-color : white; border-radius: " + self.border_curvature + "; color: rgba(5, 105, 185, 255); "
                                                                                                            "}"
                                                                                                            "QListView"
                                                                                                            "{"
                                                                                                            "background-color: white; selection-color: green"
                                                                                                            "}"
                                                                                                            "QComboBox::item:selected { background-color: rgba(255, 210, 0, 255); color: white }"
                                      )
        self.comboBox_4.addItem("")
        self.comboBox_4.addItem("")
        self.comboBox_4.addItem("")
        self.comboBox_4.addItem("")
        self.comboBox_4.addItem("")
        self.comboBox_4.addItem("")
        self.comboBox_4.addItem("")
        self.comboBox_4.addItem("")
        self.comboBox_4.addItem("")
        self.comboBox_4.addItem("")
        self.comboBox_4.addItem("")
        self.comboBox_4.addItem("")
        self.pushButton_5 = QtWidgets.QPushButton(Batch_Widget)
        self.pushButton_5.setGeometry(
            QtCore.QRect(self.left_margin, self.top_margin + self.metal_selection_section_gap + self.option_button_height + self.button_spacing, self.option_button_widht, self.option_button_height))
        self.pushButton_5.setObjectName("pushButton_5")
        self.pushButton_5.setStyleSheet("QPushButton"
                                        "{"
                                        "background-color : white; border-radius: 10px; color: rgba(5, 105, 185, 255)"
                                        "}"
                                        "QPushButton::pressed"
                                        "{"
                                        "background-color : transparent; color: white; border-style: outset; border-width: 2px; border-radius: " + self.border_curvature + "; border-color: white"
                                                                                                                                                                           "}")

        self.lineEdit = QtWidgets.QLineEdit(Batch_Widget)
        self.lineEdit.setGeometry(
            QtCore.QRect(self.left_margin + self.option_button_widht + self.button_spacing, self.top_margin + (2 * self.metal_selection_section_gap) + self.button_spacing + self.option_button_height,
                         self.option_button_widht, self.option_button_height))
        self.lineEdit.setObjectName("lineEdit")
        self.lineEdit.setText("Output Directory")
        self.lineEdit.setAlignment(QtCore.Qt.AlignHCenter)
        self.lineEdit.setStyleSheet("background-color : white; border-radius: " + self.border_curvature + "; color: rgba(5, 105, 185, 255)")

        self.lineEdit_2 = QtWidgets.QLineEdit(Batch_Widget)
        self.lineEdit_2.setGeometry(
            QtCore.QRect(self.left_margin, self.top_margin + (2 * self.metal_selection_section_gap) + self.button_spacing + self.option_button_height,
                         self.option_button_widht, self.option_button_height))
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.lineEdit_2.setText("Output Directory")
        self.lineEdit_2.setAlignment(QtCore.Qt.AlignHCenter)
        self.lineEdit_2.setStyleSheet("background-color : white; border-radius: " + self.border_curvature + "; color: rgba(5, 105, 185, 255)")

        self.frame_4 = QtWidgets.QFrame(Batch_Widget)
        self.frame_4.setGeometry(QtCore.QRect(400, 30, 371, 541))
        self.frame_4.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_4.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_4.setObjectName("frame_4")
        self.frame_4.setStyleSheet("background-color : transparent;")
        self.tableView_2 = QtWidgets.QListWidget(self.frame_4)  # QListWidget
        self.tableView_2.setGeometry(QtCore.QRect(36, round(self.viewer_height + 22.75 + 22.75), self.viewer_width, self.viewer_height))
        self.tableView_2.setObjectName("tableView_2")
        self.tableView_2.setStyleSheet(
            "background-color : transparent; border-style: outset; border-width: " + self.border_width + "; border-radius: " + self.border_curvature + "; border-color: white; color: white; ")
        self.tableWidget = QtWidgets.QTableWidget(self.frame_4)
        self.tableWidget.setGeometry(QtCore.QRect(36, 23, self.viewer_width, self.viewer_height))
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.verticalHeader().hide()
        self.tableWidget.setStyleSheet("QTableWidget"
                                       "{"
                                       "background-color: transparent; border-style: outset; border-width: " + self.border_width + "; border-radius: " + self.border_curvature + "; border-color: white; color: white;"
                                                                                                                                                                                 "}"
                                                                                                                                                                                 "QHeaderView::section"
                                                                                                                                                                                 "{"
                                                                                                                                                                                 "background-color: white; color: rgba(5, 105, 185, 255);"
                                                                                                                                                                                 "}"
                                       )

        self.tableWidget.setColumnCount(3)
        self.tableWidget.setRowCount(0)
        self.tableWidget.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.tableWidget.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(2, item)
        self.top_image_background = QtWidgets.QFrame(self.frame_4)
        self.top_image_background.setGeometry(36, round((self.viewer_height * 2) + 22.75 + 22.75 + 22.75), self.viewer_width, self.viewer_height)
        self.top_image_background.setStyleSheet(
            "background-color: black; border-radius: 20px; border-style: outset; border-width: " + self.border_width + "; border-radius: " + self.border_curvature + "; border-color: white; ")

        self.svgWidget = QtSvg.QSvgWidget('/Users/cianclarke/Documents/PhD/Complex_Assembly/GUI/topology_images/blank.svg',
                                          self.frame_4)  # The fact that I put in batch widget means that i add it to that widget
        self.svgWidget.setGeometry(-25, 375, 350, 130)

        self.frame_4.raise_()
        self.frame_2.raise_()
        self.pushButton_3.raise_()
        self.comboBox_4.raise_()
        self.pushButton_5.raise_()
        self.pushButton_4.raise_()
        self.lineEdit.raise_()
        self.frame.raise_()
        self.Save.raise_()
        self.comboBox.raise_()
        self.comboBox_2.raise_()
        self.comboBox_3.raise_()
        self.pushButton.raise_()
        self.pushButton_2.raise_()
        QtCore.QMetaObject.connectSlotsByName(Batch_Widget)
        _translate = QtCore.QCoreApplication.translate
        Batch_Widget.setWindowTitle(_translate("Batch_Widget", "Batch"))
        self.Save.setText(_translate("Batch_Widget", "Save"))
        self.comboBox_2.setItemText(0, _translate("Batch_Widget", "Metal"))
        self.comboBox_2.setItemText(1, _translate("Batch_Widget", "Sc"))
        self.comboBox_2.setItemText(2, _translate("Batch_Widget", "Ti"))
        self.comboBox_2.setItemText(3, _translate("Batch_Widget", "V"))
        self.comboBox_2.setItemText(4, _translate("Batch_Widget", "Cr"))
        self.comboBox_2.setItemText(5, _translate("Batch_Widget", "Mn"))
        self.comboBox_2.setItemText(6, _translate("Batch_Widget", "Fe"))
        self.comboBox_2.setItemText(7, _translate("Batch_Widget", "Co"))
        self.comboBox_2.setItemText(8, _translate("Batch_Widget", "Ni"))
        self.comboBox_2.setItemText(9, _translate("Batch_Widget", "Cu"))
        self.comboBox_2.setItemText(10, _translate("Batch_Widget", "Zn"))
        self.comboBox_2.setItemText(11, _translate("Batch_Widget", "Ir"))
        self.comboBox_2.setItemText(12, _translate("Batch_Widget", "Ru"))
        self.comboBox_3.setItemText(0, _translate("Batch_Widget", "Charge"))
        self.comboBox_3.setItemText(1, _translate("Batch_Widget", "0"))
        self.comboBox_3.setItemText(2, _translate("Batch_Widget", "+1"))
        self.comboBox_3.setItemText(3, _translate("Batch_Widget", "+2"))
        self.comboBox_3.setItemText(4, _translate("Batch_Widget", "+3"))
        self.comboBox_3.setItemText(5, _translate("Batch_Widget", "+4"))
        self.comboBox_3.setItemText(6, _translate("Batch_Widget", "+5"))
        self.comboBox_3.setItemText(7, _translate("Batch_Widget", "+6"))
        self.comboBox_3.setItemText(8, _translate("Batch_Widget", "+7"))
        self.comboBox.setItemText(0, _translate("Batch_Widget", "Spin"))
        self.comboBox.setItemText(1, _translate("Batch_Widget", "High"))
        self.comboBox.setItemText(2, _translate("Batch_Widget", "Low"))
        self.pushButton.setText(_translate("Batch_Widget", "Add Metal"))
        self.pushButton_2.setText(_translate("Batch_Widget", "Remove Metal"))
        self.pushButton_3.setText(_translate("Batch_Widget", "Add"))
        self.pushButton_4.setText(_translate("Batch_Widget", "Remove"))
        self.comboBox_4.setItemText(0, _translate("Batch_Widget", "Select Topology"))
        self.comboBox_4.setItemText(1, _translate("Batch_Widget", "[2, 1, 0]"))
        self.comboBox_4.setItemText(2, _translate("Batch_Widget", "[2, 1, 1, [\"1\", \"2\", \"2\"]]"))
        self.comboBox_4.setItemText(3, _translate("Batch_Widget", "[2, 1, 1, [\"1\", \"2\", \"3\"]]"))
        self.comboBox_4.setItemText(4, _translate("Batch_Widget", "[ 2, 2, [ \"1\", \"1\" ] ]"))
        self.comboBox_4.setItemText(5, _translate("Batch_Widget", "[ 2, 2, [ \"1\", \"2\" ] ]"))
        self.comboBox_4.setItemText(6, _translate("Batch_Widget", "[3, 2, 0]"))
        self.comboBox_4.setItemText(7, _translate("Batch_Widget", "[4, 1, 0]"))
        self.comboBox_4.setItemText(8, _translate("Batch_Widget", "[4, 1, 1, [\"1\", \"2\", \"2\"]]"))
        self.comboBox_4.setItemText(9, _translate("Batch_Widget", "[4, 1, 1, [\"1\", \"2\", \"3\"]]"))
        self.comboBox_4.setItemText(10, _translate("Batch_Widget", "[5, 0]"))
        self.comboBox_4.setItemText(11, _translate("Batch_Widget", "[5, 1]"))
        self.pushButton_5.setText(_translate("Batch_Widget", "View"))

        item = self.tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("Batch_Widget", "Metal"))
        item = self.tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("Batch_Widget", "Oxidation State"))
        item = self.tableWidget.horizontalHeaderItem(2)
        item.setText(_translate("Batch_Widget", "Spin"))
        self.tableWidget.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.lineEdit_3 = QtWidgets.QLineEdit(Batch_Widget)
        self.lineEdit_3.setObjectName(u"lineEdit_3")
        self.lineEdit_3.setGeometry(QtCore.QRect(self.left_margin + self.button_spacing + self.option_button_widht,
                                                 self.top_margin + (2 * self.metal_selection_section_gap) + (2 * self.button_spacing) + (2 * self.option_button_height), (self.option_button_widht),
                                                 self.option_button_height))
        self.lineEdit_2.setText("Input Directory")
        self.lineEdit_3.setText("Name Batch")
        self.lineEdit_3.setAlignment(QtCore.Qt.AlignHCenter)
        self.lineEdit_3.setStyleSheet("background-color : white; border-radius: " + self.border_curvature + "; color: rgba(5, 105, 185, 255)")
        self.radioButton = QtWidgets.QRadioButton(Batch_Widget)
        self.radioButton.setStyleSheet("background-color : black; color: white;")
        self.radioButton.setGeometry(
            QtCore.QRect(self.left_margin, self.top_margin + (2 * self.metal_selection_section_gap) + (2 * self.button_spacing) + (2 * self.option_button_height), (self.option_button_widht),
                         self.option_button_height))
        self.radioButton.setText("Optimise - UFF")
        self.radioButton.setStyleSheet("QRadioButton"
                                       "{"
                                       "background-color: transparent;"
                                       "color: white;"
                                       "}"

                                       "QRadioButton::indicator::unchecked"
                                       "{"
                                       "width: 16px;"
                                       "height: 16px;"
                                       "border: 2px solid;"
                                       "border-radius: 8px;"
                                       "border-color: black;"
                                       "background-color: rgba(0, 180, 170, 255)"

                                       "}"

                                       "QRadioButton::indicator::checked"
                                       "{"
                                       "width: 16px;"
                                       "height: 16px;"
                                       "border: 2px solid;"
                                       "border-radius: 8px;"
                                       "border-color: black;"
                                       "background-color: rgba(255, 210, 0, 255)"

                                       "}"
                                       )
        self.spinBox = QtWidgets.QSpinBox(Batch_Widget)
        self.spinBox.setObjectName(u"spinBox")
        self.spinBox.setGeometry(QtCore.QRect(self.left_margin + self.option_button_widht + self.button_spacing,
                                              self.top_margin + (2 * self.metal_selection_section_gap) + (3 * self.button_spacing) + (3 * self.option_button_height), self.option_button_widht,
                                              self.option_button_height))
        self.spinBox.setMaximum(999999999)
        self.spinBox.setStyleSheet("QSpinBox "
                                   "{"
                                   "padding-right: 20px;"
                                   "background-color : white; border-radius: " + self.border_curvature + "; color: rgba(5, 105, 185, 255)"
                                                                                                         "}"

                                                                                                         "QSpinBox::down-button"
                                                                                                         "{"
                                                                                                         "subcontrol-origin: border;"
                                                                                                         "subcontrol-position: bottom right;"
                                                                                                         "width: 12px;"
                                                                                                         "height: 12px;"
                                                                                                         "border-width: 1px;"
                                                                                                         "background-color:rgba(255, 210, 0, 255);"
                                                                                                         "left: -10px;"
                                                                                                         "border-radius: 6px;"
                                                                                                         "}"

                                                                                                         "QSpinBox::down-button:pressed"
                                                                                                         "{"
                                                                                                         "background-color:rgba(240, 184, 0, 255);"
                                                                                                         "}"
                                                                                                         "QSpinBox::up-button:pressed"
                                                                                                         "{"
                                                                                                         "background-color:rgba(240, 184, 0, 255);"
                                                                                                         "}"

                                                                                                         "QSpinBox::up-button"
                                                                                                         "{"
                                                                                                         "subcontrol-origin: border;"
                                                                                                         "subcontrol-position: top right;"
                                                                                                         "width: 12px;"
                                                                                                         "height: 12px;"
                                                                                                         "border-width: 1px;"
                                                                                                         "background-color:rgba(255, 210, 0, 255);"
                                                                                                         "left: -10px;"
                                                                                                         "border-radius: 6px;"
                                                                                                         "}"
                                   )
        self.spinBox.setSpecialValueText("MAX. Complex")

        self.spinBox_2 = QtWidgets.QSpinBox(Batch_Widget)
        self.spinBox_2.setObjectName(u"spinBox_2")
        self.spinBox_2.setGeometry(QtCore.QRect(self.left_margin, self.top_margin + self.metal_selection_section_gap, self.option_button_widht,
                                                self.option_button_height))
        self.spinBox_2.setMaximum(999999999)
        self.spinBox_2.setStyleSheet("QSpinBox "
                                     "{"
                                     "padding-right: 20px;"
                                     "background-color : white; border-radius: " + self.border_curvature + "; color: rgba(5, 105, 185, 255)"
                                                                                                           "}"

                                                                                                           "QSpinBox::down-button"
                                                                                                           "{"
                                                                                                           "subcontrol-origin: border;"
                                                                                                           "subcontrol-position: bottom right;"
                                                                                                           "width: 12px;"
                                                                                                           "height: 12px;"
                                                                                                           "border-width: 1px;"
                                                                                                           "background-color:rgba(255, 210, 0, 255);"
                                                                                                           "left: -10px;"
                                                                                                           "border-radius: 6px;"
                                                                                                           "}"

                                                                                                           "QSpinBox::down-button:pressed"
                                                                                                           "{"
                                                                                                           "background-color:rgba(240, 184, 0, 255);"
                                                                                                           "}"
                                                                                                           "QSpinBox::up-button:pressed"
                                                                                                           "{"
                                                                                                           "background-color:rgba(240, 184, 0, 255);"
                                                                                                           "}"

                                                                                                           "QSpinBox::up-button"
                                                                                                           "{"
                                                                                                           "subcontrol-origin: border;"
                                                                                                           "subcontrol-position: top right;"
                                                                                                           "width: 12px;"
                                                                                                           "height: 12px;"
                                                                                                           "border-width: 1px;"
                                                                                                           "background-color:rgba(255, 210, 0, 255);"
                                                                                                           "left: -10px;"
                                                                                                           "border-radius: 6px;"
                                                                                                           "}"
                                     )
        self.spinBox_2.setSpecialValueText("Random Seed")

        self.label_10 = QtWidgets.QLabel(Batch_Widget)
        self.label_10.setObjectName(u"label_10")
        self.label_10.setGeometry(
            QtCore.QRect(self.left_margin, self.top_margin + (2 * self.metal_selection_section_gap) + (3 * self.button_spacing) + (3 * self.option_button_height), self.option_button_widht,
                         self.option_button_height))
        self.label_10.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.label_10.setStyleSheet("background-color : transparent; color: white; border-style: outset; border-width: 2px; border-radius: " + self.border_curvature + "; border-color: white")
        self.label_10.setText(QtCore.QCoreApplication.translate("Batch_Widget", u"Max.  Complex", None))

        self.comboBox_isomer = QtWidgets.QComboBox(Batch_Widget)
        self.comboBox_isomer.setObjectName(u"combo_box_isomer")
        self.comboBox_isomer.setGeometry(
            QtCore.QRect(self.left_margin, self.top_margin + self.metal_selection_section_gap + (2 * self.option_button_height) + (2 * self.button_spacing), self.option_button_widht,
                         self.option_button_height))
        self.comboBox_isomer.addItem("Isomers")
        self.comboBox_isomer.addItem("Generate All")
        self.comboBox_isomer.addItem("Generate Lowest Energy")
        self.comboBox_isomer.setStyleSheet("QComboBox::drop-down { "
                                           "subcontrol-origin: padding; "
                                           "subcontrol-position: top right; "
                                           "width: 10px; "
                                           "border-left-width: 5px; "
                                           "border-left-color: rgba(255, 210, 0, 255); "
                                           "border-left-style: solid; "
                                           "border-top-right-radius: 3px; "
                                           "border-bottom-right-radius: 3px; "
                                           "}"
                                           "QComboBox"
                                           "{"
                                           "background-color : white; border-radius: " + self.border_curvature + "; color: rgba(5, 105, 185, 255); "
                                                                                                                 "}"
                                                                                                                 "QListView"
                                                                                                                 "{"
                                                                                                                 "background-color: white; selection-color: green"
                                                                                                                 "}"
                                                                                                                 "QComboBox::item:selected { background-color: rgba(255, 210, 0, 255); color: white }"
                                           )
        ###clicked buttons###                                                       #
        self.pushButton_3.clicked.connect(self.write_topologies_to_dialog_box)  #
        self.pushButton_4.clicked.connect(self.remove_topologies_from_dialog_box)  #
        self.pushButton.clicked.connect(self.write_metals_to_dialog_box)  #
        self.pushButton_2.clicked.connect(self.remove_metals_from_dialog_box)  #
        self.pushButton_5.clicked.connect(self.view_topologies)
        #####################                                                       #

    def close_assembly_widget(self):  # This function closes the widget
        if (self.lineEdit_3.text() != "Name Batch") and ((self.lineEdit.text() != "Output Directory") and (self.lineEdit.text() != "")) and (
                self.spinBox.text() != "MAX. Complex") and self.tableView_2.count() and self.tableWidget.rowCount() and (str(self.comboBox_isomer.currentText()) != "Isomers") and (
                (self.lineEdit_2.text() != "Input Directory") and (self.lineEdit_2.text() != "")) and (self.spinBox_2.text() != "Random Seed"):
            self.tmp_Dic.clear()
            self.tmp_Dic.update({"Name": str(self.lineEdit_3.text())})
            self.tmp_Dic.update({"Output_Path": self.lineEdit.text()})
            self.tmp_Dic.update({"Input_Path": self.lineEdit_2.text()})

            self.tmp_Dic.update({"Input_Path": self.lineEdit_2.text()})

            self.tmp_Dic.update({"MAX_num_complexes": self.spinBox.text()})
            self.tmp_Dic.update({"Isomers": self.comboBox_isomer.currentText()})
            self.tmp_Dic.update({"Optimisation_Choice": str(self.radioButton.isChecked())})
            self.tmp_Dic.update({"Random_Seed": str(self.spinBox_2.text())})
            for i in range(self.tableView_2.count()):
                self.tmp_Dic.update({"Topology_{}".format(str(i + 1)): self.tableView_2.item(i).text()})

            for j in range(self.tableWidget.rowCount()):
                self.tmp_Dic.update({"Metal_{}".format(str(j + 1)): [self.tableWidget.item(j, 0).text(), self.tableWidget.item(j, 1).text(), self.tableWidget.item(j, 2).text()]})
            assembly_input_list.append(self.tmp_Dic)
            print(assembly_input_list)
            Window.update_batch_combo_box(ui)  # ui is the instance of the other class we are working with
            self.close()
        else:
            self.Save.setStyleSheet("QPushButton"
                                    "{"
                                    "background-color : orange; border-radius: 15px; color: black; font-size: 15pt"
                                    "}")
            self.Save.setText("!!!WARNING COMPLETE ALL FIELDS!!!")

        print(json.dumps(self.tmp_Dic, indent=4))

    def mousePressEvent(self, a0: QtGui.QMouseEvent) -> None:
        self.Save.setStyleSheet("")
        self.Save.setText("Save")
        self.Save.setStyleSheet("QPushButton"
                                "{"
                                "background-color : rgba(255, 210, 0, 255); border-radius: 15px; color: rgba(5, 105, 185, 255); font-size: 20pt"
                                "}"
                                "QPushButton::pressed"
                                "{"
                                "background-color : transparent; color: rgba(255, 210, 0, 255); border-style: outset; border-width: 2px; border-radius: 15px; border-color: rgba(255, 210, 0, 255)"
                                "}")

    def write_topologies_to_dialog_box(self):
        topology_choice = self.comboBox_4.currentText()
        if topology_choice != "Select Topology":
            self.tableView_2.addItem(str(topology_choice))

    def remove_topologies_from_dialog_box(self):
        value = self.tableView_2.currentRow()
        self.tableView_2.takeItem(value)

    def iterAllItems(self):  # This method will extract all items from a qlistwidget
        self.all_items = self.tableView_2.findItems('', QtCore.Qt.MatchFlag)
        print(self.all_items)

    def write_metals_to_dialog_box(self):
        metal_choice = self.comboBox_2.currentText()
        oxidation_state_choice = self.comboBox_3.currentText()
        spin_choice = self.comboBox.currentText()
        if metal_choice != "Metal" and oxidation_state_choice != "Charge" and spin_choice != "Spin":
            self.current_num = self.tableWidget.rowCount()
            print(self.current_num + 1)
            self.tableWidget.setRowCount(self.current_num + 1)
            self.tableWidget.setItem(self.current_num, 0, QTableWidgetItem(str(metal_choice)))
            self.tableWidget.setItem(self.current_num, 1, QTableWidgetItem(str(oxidation_state_choice)))
            self.tableWidget.setItem(self.current_num, 2, QTableWidgetItem(str(spin_choice)))

    def remove_metals_from_dialog_box(self):
        self.current_row = self.tableWidget.currentRow()
        self.tableWidget.removeRow(self.current_row)

    def view_topologies(self):
        self.current_topology = self.comboBox_4.currentText()
        self.svgWidget.close()

        print(self.current_topology)
        if self.current_topology == "[2, 1, 0]":

            self.svgWidget = QtSvg.QSvgWidget('/Users/cianclarke/Documents/PhD/Complex_Assembly/GUI/topology_images/2_1_1__1_2_3.svg',
                                              self.frame_4)  # The fact that I put in batch widget means that i add it to that widget
            self.svgWidget.setGeometry(115, 380, 130, 130)
            self.svgWidget.show()

        elif self.current_topology == "[2, 1, 1, [\"1\", \"2\", \"2\"]]":
            self.svgWidget = QtSvg.QSvgWidget('/Users/cianclarke/Documents/PhD/Complex_Assembly/GUI/topology_images/2_1_1__1_2_2.svg',
                                              self.frame_4)  # The fact that I put in batch widget means that i add it to that widget
            self.svgWidget.setGeometry(115, 380, 130, 130)
            self.svgWidget.show()


        elif self.current_topology == "[2, 1, 1, [\"1\", \"2\", \"3\"]]":
            self.svgWidget = QtSvg.QSvgWidget('/Users/cianclarke/Documents/PhD/Complex_Assembly/GUI/topology_images/2_1_1__1_2_3.svg',
                                              self.frame_4)  # The fact that I put in batch widget means that i add it to that widget
            self.svgWidget.setGeometry(115, 380, 130, 130)
            self.svgWidget.show()


        elif self.current_topology == "[ 2, 2, [ \"1\", \"1\" ] ]":
            self.svgWidget = QtSvg.QSvgWidget('/Users/cianclarke/Documents/PhD/Complex_Assembly/GUI/topology_images/2_2__1_1.svg',
                                              self.frame_4)  # The fact that I put in batch widget means that i add it to that widget
            self.svgWidget.setGeometry(115, 380, 130, 130)
            self.svgWidget.show()


        elif self.current_topology == "[ 2, 2, [ \"1\", \"2\" ] ]":
            self.svgWidget = QtSvg.QSvgWidget('/Users/cianclarke/Documents/PhD/Complex_Assembly/GUI/topology_images/2_2__1_2.svg',
                                              self.frame_4)  # The fact that I put in batch widget means that i add it to that widget
            self.svgWidget.setGeometry(115, 380, 130, 130)
            self.svgWidget.show()


        elif self.current_topology == "[3, 2, 0]":
            self.svgWidget = QtSvg.QSvgWidget('/Users/cianclarke/Documents/PhD/Complex_Assembly/GUI/topology_images/3_2_0.svg',
                                              self.frame_4)  # The fact that I put in batch widget means that i add it to that widget
            self.svgWidget.setGeometry(115, 380, 130, 130)
            self.svgWidget.show()


        elif self.current_topology == "[4, 1, 0]":
            self.svgWidget = QtSvg.QSvgWidget('/Users/cianclarke/Documents/PhD/Complex_Assembly/GUI/topology_images/4_1_1__1_2_3.svg',
                                              self.frame_4)  # The fact that I put in batch widget means that i add it to that widget
            self.svgWidget.setGeometry(115, 380, 130, 130)
            self.svgWidget.show()


        elif self.current_topology == "[4, 1, 1, [\"1\", \"2\", \"2\"]]":
            self.svgWidget = QtSvg.QSvgWidget('/Users/cianclarke/Documents/PhD/Complex_Assembly/GUI/topology_images/4_1_1__1_2_2.svg',
                                              self.frame_4)  # The fact that I put in batch widget means that i add it to that widget
            self.svgWidget.setGeometry(115, 380, 130, 130)
            self.svgWidget.show()


        elif self.current_topology == "[4, 1, 1, [\"1\", \"2\", \"3\"]]":
            self.svgWidget = QtSvg.QSvgWidget('/Users/cianclarke/Documents/PhD/Complex_Assembly/GUI/topology_images/4_1_1__1_2_3.svg',
                                              self.frame_4)  # The fact that I put in batch widget means that i add it to that widget
            self.svgWidget.setGeometry(115, 380, 130, 130)
            self.svgWidget.show()


        elif self.current_topology == "[5, 0]":
            self.svgWidget = QtSvg.QSvgWidget('/Users/cianclarke/Documents/PhD/Complex_Assembly/GUI/topology_images/5_1.svg',
                                              self.frame_4)  # The fact that I put in batch widget means that i add it to that widget
            self.svgWidget.setGeometry(115, 380, 130, 130)
            self.svgWidget.show()


        elif self.current_topology == "[5, 1]":
            self.svgWidget = QtSvg.QSvgWidget('/Users/cianclarke/Documents/PhD/Complex_Assembly/GUI/topology_images/5_1.svg',
                                              self.frame_4)  # The fact that I put in batch widget means that i add it to that widget
            self.svgWidget.setGeometry(115, 380, 130, 130)
            self.svgWidget.show()

        else:
            self.svgWidget.close()


class Window(QtWidgets.QMainWindow):

    def main_window(self, MainWindow):
        self.setFixedSize(800, 600)
        self.thread = {}  # can add threads as its values

        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.centralwidget.setStyleSheet("background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 rgba(5, 105, 185, 255), stop:1 rgba(0, 170, 205, 255));")

        self.DialogBox = QtWidgets.QListWidget(self.centralwidget)
        self.DialogBox.setGeometry(QtCore.QRect(380, 20, 391, 500))
        self.DialogBox.setObjectName("DialogBox")
        self.DialogBox.setStyleSheet("background-color: black; border-radius: 20px; color: white")
        self.DialogBox.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        self.NewBatch = QtWidgets.QPushButton(self.centralwidget)
        self.NewBatch.setGeometry(QtCore.QRect(20, 20, 131, 50))
        self.NewBatch.setObjectName("NewBatch")
        self.NewBatch.clicked.connect(self.window2)
        self.NewBatch.setStyleSheet("QPushButton"
                                    "{"
                                    "background-color : white; border-radius: 20px; color: rgba(5, 105, 185, 255); font-size: 20pt"
                                    "}"
                                    "QPushButton::pressed"
                                    "{"
                                    "background-color : transparent; color: white; border-style: outset; border-width: 2px; border-radius: 20px; border-color: white"
                                    "}")

        self.Run = QtWidgets.QPushButton(self.centralwidget)
        self.Run.setGeometry(QtCore.QRect(20, 480, 131, 50))
        self.Run.setObjectName("Run")
        self.Run.clicked.connect(self.Run_Program)
        self.Run.setStyleSheet("QPushButton"
                               "{"
                               "background-color : rgba(255, 210, 0, 255); border-radius: 20px; color: rgba(5, 105, 185, 255); font-size: 30pt"
                               "}"
                               "QPushButton::pressed"
                               "{"
                               "background-color : transparent; color: rgba(255, 210, 0, 255); border-style: outset; border-width: 2px; border-radius: 20px; border-color: rgba(255, 210, 0, 255)"
                               "}")

        self.Cancel = QtWidgets.QPushButton(self.centralwidget)
        self.Cancel.setGeometry(QtCore.QRect(190, 480, 131, 50))
        self.Cancel.setObjectName("Cancel")
        self.Cancel.clicked.connect(self.close_GUI)
        self.Cancel.setStyleSheet("QPushButton"
                                  "{"
                                  "background-color : rgba(0, 180, 170, 255); border-radius: 20px; color: rgba(5, 105, 185, 255); font-size: 30pt"
                                  "}"
                                  "QPushButton::pressed"
                                  "{"
                                  "background-color : transparent; color: rgba(0, 180, 170, 255); border-style: outset; border-width: 2px; border-radius: 20px; border-color: rgba(0, 180, 170, 255)"
                                  "}")

        self.ViewBatch = QtWidgets.QPushButton(self.centralwidget)
        self.ViewBatch.setGeometry(QtCore.QRect(20, 90, 131, 50))
        self.ViewBatch.setObjectName("ViewBatch")
        self.ViewBatch.clicked.connect(self.view_batch)
        self.ViewBatch.setStyleSheet("QPushButton"
                                     "{"
                                     "background-color : white; border-radius: 20px; color: rgba(5, 105, 185, 255); font-size: 20pt"
                                     "}"
                                     "QPushButton::pressed"
                                     "{"
                                     "background-color : transparent; color: white; border-style: outset; border-width: 2px; border-radius: 20px; border-color: white"
                                     "}")

        self.DeleteBatch = QtWidgets.QPushButton(self.centralwidget)
        self.DeleteBatch.setGeometry(QtCore.QRect(20, 160, 131, 50))
        self.DeleteBatch.setObjectName("DeleteBatch")
        self.DeleteBatch.clicked.connect(self.delete_batch)
        self.DeleteBatch.setStyleSheet("QPushButton"
                                       "{"
                                       "background-color : white; border-radius: 20px; color: rgba(5, 105, 185, 255); font-size: 20pt"
                                       "}"
                                       "QPushButton::pressed"
                                       "{"
                                       "background-color : transparent; color: white; border-style: outset; border-width: 2px; border-radius: 20px; border-color: white"
                                       "}")

        self.VersionLabel = QtWidgets.QLabel(self.centralwidget)
        self.VersionLabel.setGeometry(QtCore.QRect(380, 540, 391, 16))
        self.VersionLabel.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.VersionLabel.setObjectName("VersionLabel")
        self.VersionLabel.setStyleSheet("background-color : transparent; color: white; font-size: 15pt")

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        self.menubar.setObjectName("menubar")

        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.NewBatch.setText(_translate("MainWindow", "New"))
        self.Run.setText(_translate("MainWindow", "O"))
        self.Cancel.setText(_translate("MainWindow", "X"))
        self.ViewBatch.setText(_translate("MainWindow", "View"))
        self.VersionLabel.setText(_translate("MainWindow", "V 1.0.0"))
        self.DeleteBatch.setText(_translate("MainWindow", "Delete"))

        self.comboBox = QtWidgets.QComboBox(self.centralwidget)
        self.comboBox.setObjectName(u"batch_combo_box")
        self.comboBox.setGeometry(QtCore.QRect(190, 20, 131, 51))
        self.comboBox.addItem("   Select")
        self.comboBox.setStyleSheet("QComboBox::drop-down { "
                                    "subcontrol-origin: padding; "
                                    "subcontrol-position: top right; "
                                    "width: 22px; "
                                    "border-left-width: 5px; "
                                    "border-left-color: rgba(255, 210, 0, 255); "
                                    "border-left-style: solid; "
                                    "border-top-right-radius: 3px; "
                                    "border-bottom-right-radius: 3px; "
                                    "}"
                                    "QComboBox"
                                    "{"
                                    "background-color : white; border-radius: 20px; color: rgba(5, 105, 185, 255); font-size: 20pt"
                                    "}"
                                    "QListView"
                                    "{"
                                    "background-color: white; selection-color: green"
                                    "}"
                                    "QComboBox::item:selected { background-color: rgba(255, 210, 0, 255); color: white }"
                                    )

    def window2(self):  # <===
        self.w = Window2()
        self.w.setupUi(self.w)
        self.w.show()

    def view_batch(self):
        if self.comboBox.currentText() == "   Select":
            pass
        else:
            print(self.comboBox.currentText())
            for batch in assembly_input_list:
                if batch["Name"] == str(self.comboBox.currentText()).replace("   ", ""):
                    print("im here")
                    for key, value in batch.items():
                        self.DialogBox.addItem("   " + str(key) + ": " + str(value))
                else:
                    pass
            self.DialogBox.addItem(" ")
            self.DialogBox.addItem("   ############")
            self.DialogBox.addItem(" ")

    def update_batch_combo_box(self):
        print("clicked")
        self.comboBox.addItem("   " + str(assembly_input_list[-1]["Name"]))

    def delete_batch(self):
        deletion_index = "var"
        selected_batch = str(self.comboBox.currentText()).replace("   ", "")
        i = 0
        for batch in assembly_input_list:
            if str(batch["Name"]) == str(selected_batch):
                deletion_index = i
            else:
                pass
            i = i + 1
        if (deletion_index == "var") and (self.comboBox.currentText() != "   Select"):
            print("!!!Fatal Batch Deletion Error!!!")
            exit()
        elif (deletion_index == "var") and (self.comboBox.currentText() == "   Select"):
            pass
        else:
            assembly_input_list.pop(deletion_index)
            current_index = self.comboBox.findText(self.comboBox.currentText())
            self.comboBox.removeItem(current_index)

    def close_GUI(self):
        print("cancelling GUI")
        MainWindow.close()

    def Run_Program(self):
        self.thread[1] = ThreadClass(parent=None, index=1)
        self.thread[1].start()
        self.thread[1].any_signal.connect(self.my_function)
        self.Run.setEnabled(False)
        self.Run.setStyleSheet("QPushButton"
                               "{"
                               "background-color : transparent; color: rgba(255, 210, 0, 255); border-style: outset; border-width: 2px; border-radius: 20px; border-color: rgba(255, 210, 0, 255); font-size: 30pt"
                               "}")

    def my_function(self, counter):
        cnt = counter
        index = self.sender().index
        if index == 1:
            self.progressBar.setValue(cnt)
        if index == 2:
            self.progressBar_2.setValue(cnt)
        if index == 3:
            self.progressBar_3.setValue(cnt)


class ThreadClass(QtCore.QThread):
    any_signal = QtCore.pyqtSignal(int)

    def __init__(self, parent=None, index=0):
        super(ThreadClass, self).__init__(parent)
        self.index = index
        self.is_running = True

    def run(self):
        print('Starting thread...', self.index)
        Assembly_instance = Assembly(list_of_batches=assembly_input_list)
        Assembly_instance.assembly_main()
        self.is_running = False
        self.terminate()

    def stop(self):
        self.is_running = False
        print('Stopping thread...', self.index)
        self.terminate()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("Fusion")
    MainWindow = QtWidgets.QMainWindow()
    ui = Window()
    ui.main_window(MainWindow)
    MainWindow.show()
    sys.exit(app.exec())
