# In this app I am going to add a button that will generate xyz files for all the complexes that are coloured green

import ase
import numpy as np
import ase.visualize
import sys
from PyQt5.QtWidgets import QWidget, QSizePolicy
from PyQt5.QtGui import QColor
import pyqtgraph as pg
import json
import ast
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import QSize
import qtawesome as Qta
import os

# Note!!! this script requires packages that are not part of the DART installation so this file should not be made public
f = open('CASINI_OUTPUT_130923.json')
# f = open('/Users/cianclarke/Documents/PhD/Complex_Assembly/CreateTMC/Gaussian/CASINI_DATA/CASINI_OUTPUT.json')
NBO_dict = json.load(f)


class Window(QtWidgets.QMainWindow):
    resized = QtCore.pyqtSignal()

    def __init__(self, parent=None, input_data=None):
        super(Window, self).__init__(parent=parent)
        # Here we assign some variables to none so we don't get errors
        self.centralwidget = None
        self.pushButton = None
        self.pushButton_2a = None
        self.pushButton_2b = None
        self.pushButton_3 = None
        self.DialogBox = None
        self.graphicsView = None
        self.proxy = None
        self.MOUSE_X = None
        self.MOUSE_Y = None
        self.crosshair_v = None
        self.crosshair_h = None
        self.chosen_index = None
        self.plot = None
        self.bound_x = None
        self.bound_y = None
        self.current_plot_num = None
        self.downloadable_xyz_indexes = None  # Here we save the indexes we want to let the user download

        self.colours = []
        self.x_axis_list = []
        self.y_axis_list = []
        self.other_data = []
        self.Name = []
        self.structure = []
        self.symbols = []
        self.symbols_size = []

        # Here we load the data
        self.data = input_data

        # These are some of the parameters that control the layout of the main window
        self.window_x = 1920
        self.window_y = 1080
        self.spacing = 30
        self.widget_ratio = 0.5
        self.graph_x_length = int((self.window_x - (1.5 * self.spacing)) * self.widget_ratio)
        self.button_height = int(3 * self.spacing)

        self.primary_colour_HEX = "#313131"
        self.primary_colour_R = self.hex_to_rgb(self.primary_colour_HEX)[0]
        self.primary_colour_G = self.hex_to_rgb(self.primary_colour_HEX)[1]
        self.primary_colour_B = self.hex_to_rgb(self.primary_colour_HEX)[2]
        self.primary_colour_A = 255
        self.primary_colour = f"rgba({self.primary_colour_R},{self.primary_colour_G}, {self.primary_colour_B}, {self.primary_colour_A})"

        self.secondary_colour_HEX = "#e86f0c"
        self.secondary_colour_R = self.hex_to_rgb(self.secondary_colour_HEX)[0]
        self.secondary_colour_G = self.hex_to_rgb(self.secondary_colour_HEX)[1]
        self.secondary_colour_B = self.hex_to_rgb(self.secondary_colour_HEX)[2]
        self.secondary_colour_A = 255
        self.secondary_colour = f"rgba({self.secondary_colour_R},{self.secondary_colour_G}, {self.secondary_colour_B}, {self.secondary_colour_A})"

        self.tertiary_colour_HEX = "#62a68b"
        self.tertiary_colour_R = self.hex_to_rgb(self.tertiary_colour_HEX)[0]
        self.tertiary_colour_G = self.hex_to_rgb(self.tertiary_colour_HEX)[1]
        self.tertiary_colour_B = self.hex_to_rgb(self.tertiary_colour_HEX)[2]
        self.tertiary_colour_A = 255
        self.tertiary_colour = f"rgba({self.tertiary_colour_R},{self.tertiary_colour_G}, {self.tertiary_colour_B}, {self.tertiary_colour_A})"

        # set up the ui
        self.setupUi(self)

        # Initiate the plot
        start_ = 0
        self.current_plot_number = start_
        self.plot_list = [self.initiate_plot]
        self.plot_list[start_]()

        # if the window resizes we execute some function
        self.resized.connect(self.someFunction)

    def setupUi(self, MainWindow):

        # Main Window Parameters
        MainWindow.setObjectName("MainWindow")
        MainWindow.setWindowTitle("High-Throughput Screening Analyzer")
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(self.window_x, self.window_y)

        # Set central widget this has the gradient background colour
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")  # 64 100
        self.centralwidget.setStyleSheet(f"background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 {self.primary_colour}, stop:1 {self.primary_colour});")
        self.setMinimumSize(720, 480)
        MainWindow.setCentralWidget(self.centralwidget)

        # Button to Auto Range the graph
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setObjectName("pushButton")
        self.pushButton.clicked.connect(self.resetView)
        self.pushButton.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.pushButton.setText("AutoRange")
        self.pushButton.setStyleSheet("QPushButton"
                                      "{"
                                      f"background-color : {self.secondary_colour}; border-radius: 20px; color: white; font-size: 30pt"
                                      "}"
                                      "QPushButton::pressed"
                                      "{"
                                      f"background-color : transparent; color: {self.secondary_colour}; border-style: outset; border-width: 2px; border-radius: 20px; border-color: {self.secondary_colour}"
                                      "}")

        styling_icon = Qta.icon('fa5s.music',
                                active='fa5s.balance-scale',
                                color='blue',
                                color_active='orange')

        icon_graph = Qta.icon("msc.graph-line", color="#FFFFFF")
        self.pushButton.setIcon(icon_graph)
        self.pushButton.setIconSize(QSize(50, 50))

        # Button to cycle back
        self.pushButton_2a = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_2a.setObjectName("pushButton_2")
        self.pushButton_2a.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # todo: assign clear function somewhere else
        self.pushButton_2a.clicked.connect(self.decrement)
        self.pushButton_2a.setText("")
        self.pushButton_2a.setStyleSheet("QPushButton"
                                         "{"
                                         f"background-color : {self.secondary_colour}; border-radius: 20px; color: white; font-size: 30pt"
                                         "}"
                                         "QPushButton::pressed"
                                         "{"
                                         f"background-color : transparent; color: {self.secondary_colour}; border-style: outset; border-width: 2px; border-radius: 20px; border-color: {self.secondary_colour}"
                                         "}")

        icon_not_used = Qta.icon("mdi.arrow-left-drop-circle", color="#FFFFFF")
        self.pushButton_2a.setIcon(icon_not_used)
        self.pushButton_2a.setIconSize(QSize(50, 50))

        # Button to cycle foward
        self.pushButton_2b = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_2b.setObjectName("pushButton_2")
        self.pushButton_2b.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.pushButton_2b.clicked.connect(self.increment)
        self.pushButton_2b.setText("")
        self.pushButton_2b.setStyleSheet("QPushButton"
                                         "{"
                                         f"background-color : {self.secondary_colour}; border-radius: 20px; color: white; font-size: 30pt"
                                         "}"
                                         "QPushButton::pressed"
                                         "{"
                                         f"background-color : transparent; color: {self.secondary_colour}; border-style: outset; border-width: 2px; border-radius: 20px; border-color: {self.secondary_colour}"
                                         "}")

        icon_not_used = Qta.icon("mdi.arrow-right-drop-circle", color="#FFFFFF")
        self.pushButton_2b.setIcon(icon_not_used)
        self.pushButton_2b.setIconSize(QSize(50, 50))

        #####
        self.pushButton_2c = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_2c.setObjectName("pushButton_2")
        self.pushButton_2c.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.pushButton_2c.clicked.connect(self.view_chem_draw)
        self.pushButton_2c.setText(".XYZ")
        self.pushButton_2c.setStyleSheet("QPushButton"
                                         "{"
                                         f"background-color : {self.secondary_colour}; border-radius: 20px; color: white; font-size: 30pt"
                                         "}"
                                         "QPushButton::pressed"
                                         "{"
                                         f"background-color : transparent; color: {self.secondary_colour}; border-style: outset; border-width: 2px; border-radius: 20px; border-color: {self.secondary_colour}"
                                         "}")

        icon_not_used = Qta.icon("ph.file-arrow-down-light", color="#FFFFFF")
        self.pushButton_2c.setIcon(icon_not_used)
        self.pushButton_2c.setIconSize(QSize(50, 50))

        # Button to launch structure
        self.pushButton_3 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_3.setObjectName("pushButton_2")
        self.pushButton_3.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.pushButton_3.clicked.connect(self.Launch_Structure)
        self.pushButton_3.setText("View Structure")
        self.pushButton_3.setStyleSheet("QPushButton"
                                        "{"
                                        f"background-color : {self.secondary_colour}; border-radius: 20px; color: white; font-size: 30pt"
                                        "}"
                                        "QPushButton::pressed"
                                        "{"
                                        f"background-color : transparent; color: {self.secondary_colour}; border-style: outset; border-width: 2px; border-radius: 20px; border-color: {self.secondary_colour}"
                                        "}")

        icon_view = Qta.icon("fa5s.binoculars", color="#FFFFFF")
        self.pushButton_3.setIcon(icon_view)
        self.pushButton_3.setIconSize(QSize(50, 50))

        # This is where information is passed on to the user
        self.DialogBox = QtWidgets.QListWidget(self.centralwidget)
        self.DialogBox.setObjectName("DialogBox")
        self.DialogBox.setStyleSheet("background-color: black; border-radius: 20px; color: white")
        self.DialogBox.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        # This is where the data is plotted
        self.graphicsView = pg.PlotWidget(self.centralwidget)
        self.graphicsView.setObjectName("graphicsView")
        self.graphicsView.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.proxy = pg.SignalProxy(self.graphicsView.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMovedEvent)
        self.graphicsView.scene().sigMouseClicked.connect(self.mouseClickedEvent)
        self.add_cross_hair()

        # Connects the slots
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def initiate_plot(self):
        # Reset all the data storage
        self.current_plot_num = 1
        self.colours = []
        self.x_axis_list = []
        self.y_axis_list = []
        self.other_data = []
        self.Name = []
        self.structure = []
        self.symbols = []
        self.symbols_size = []
        self.downloadable_xyz_indexes = []

        # Plot description
        self.DialogBox.addItem(" ")
        start_banner = QtWidgets.QListWidgetItem(f"############################---[Plot:{self.current_plot_num}]---############################")
        start_banner.setForeground(QColor(self.secondary_colour_R, self.secondary_colour_G, self.secondary_colour_B, alpha=self.secondary_colour_A))
        self.DialogBox.addItem(start_banner)
        self.DialogBox.addItem("Scatter Plot of the NBO data for the Au-C and Au-N interactions")
        self.DialogBox.addItem("The white cross hair represents the Benchmark candidate complex")
        self.DialogBox.addItem("Each blue circle represents a transition metal complex from the calculations")
        end_banner = QtWidgets.QListWidgetItem(f"##################################################################")
        end_banner.setForeground(QColor(self.secondary_colour_R, self.secondary_colour_G, self.secondary_colour_B, alpha=self.secondary_colour_A))
        self.DialogBox.addItem(end_banner)
        self.DialogBox.addItem(" ")
        self.DialogBox.scrollToBottom()

        num_of_6_or_greater_metalocycles_within_bottom_left_quadrant = 0
        # Here we loop through all the complexes and decide which examples we would like to display
        for key, value in self.data.items():
            # todo !!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!
            # please be aware that this if statement is not display all the data points
            if (float(value['HTS_Data']["X-value"]) < 1000) and (float(value['HTS_Data']["Y-value"]) < 1000) and (float(value['HTS_Data']["X-value"]) > 20):
                structure_name = key
                self.x_axis_list.append(-1*float(value['HTS_Data']["X-value"]))
                self.y_axis_list.append(-1*float(value['HTS_Data']["Y-value"]))
                self.Name.append(structure_name)
                self.other_data.append(json.dumps(value['Other_Data'], indent=4))
                self.structure.append(str(value['Structure']))

                # blue --> all data points
                opacity = 0.25
                self.colours.append(pg.mkBrush(37, 188, 247, (int(opacity * 255))))
                self.symbols.append('o')
                self.symbols_size.append(10)
                if (float(value['HTS_Data']["X-value"]) < 97.77) and (float(value['HTS_Data']["Y-value"]) < 73.36):
                    num_of_6_or_greater_metalocycles_within_bottom_left_quadrant += 1

                self.downloadable_xyz_indexes.append(key)  # Here we save all the keys (indexes of relevant complexes) that we wish to download at some point

            else:
                pass

        print(num_of_6_or_greater_metalocycles_within_bottom_left_quadrant)

        # In this section of the code we need to add all of Casini's benchmark points
        self.plot_specific_data_points()

        assert len(self.x_axis_list) == len(self.y_axis_list) == len(self.colours) == len(self.symbols) == len(self.symbols_size) == len(self.Name) == len(self.other_data)
        self.graphicsView.plot(self.x_axis_list, self.y_axis_list, pen=None, symbol=self.symbols, symbolPen=None, symbolBrush=self.colours, symbolSize=self.symbols_size, )  # .setOpacity(0.25)
        self.graphicsView.setLabel(f'left', f"<span style=\"color:{self.secondary_colour_HEX};font-size:20px\">{'Au-C NBO'}</span>")
        self.graphicsView.setLabel(f'bottom', f"<span style=\"color:{self.secondary_colour_HEX};font-size:20px\">{'Au-N NBO'}</span>")
        self.bound_x = pg.InfiniteLine(angle=90, movable=False, pen=pg.mkPen(color='w', width=1.0))
        self.bound_y = pg.InfiniteLine(angle=0, movable=False, pen=pg.mkPen(color='w', width=1.0))
        self.bound_x.setPos(-1*97.77)
        self.bound_y.setPos(-1*73.36)
        #self.graphicsView.addItem(self.bound_x)
        #self.graphicsView.addItem(self.bound_y)

    def plot_specific_data_points(self):
        """
        Here we plot specific data points that we want to display
        """
        opacity_benchmark = 0.8
        r, g, b = 245, 113, 199
        symbol_size = 30

        mol_dict = [{"name": "Au-1",
                     "x": -99.11,
                     "y": -86.28,
                     "other_data": "CH2 Bridge. Original from Casini."},
                    {"name": "Au-2",
                     "x": -97.11,
                     "y": -79.67,
                     "other_data": "NH Bridge. Original from Casini."},
                    {"name": "Au-3",
                     "x": -97.77,
                     "y": -73.36,
                     "other_data": "CO Bridge. Original from Casini."},
                    {"name": "Au-4",
                     "x": -96.14,
                     "y": -57.11,
                     "other_data": "Au-4"},
                    {"name": "Au-5",
                     "x": -79.56,
                     "y": -105.19,
                     "other_data": "Au-5"},
                    {"name": "Au-6",
                     "x": -108.75,
                     "y": -49.94,
                     "other_data": "Au-6"},
                    {"name": "Au-7",
                     "x": -97.64,
                     "y": -81.34,
                     "other_data": "Au-7"},
                    ]


        for mol in mol_dict:
            # Molecule1a
            self.Name.append(mol["name"])
            self.colours.append(pg.mkBrush(r, g, b, (int(opacity_benchmark * 255))))
            self.x_axis_list.append(mol["x"])
            self.y_axis_list.append(mol["y"])
            self.other_data.append(mol["other_data"])
            self.structure.append(None)
            self.symbols.append('star')
            self.symbols_size.append(symbol_size)


    def resetView(self):
        self.graphicsView.getPlotItem().vb.enableAutoRange()

    def Launch_Structure(self):
        if self.chosen_index is not None:
            if self.structure[self.chosen_index] is not None:

                if self.chosen_index is None:
                    warning = QtWidgets.QListWidgetItem(f"!!!warning!!! -> Please select a point before launching the structure -> Proceeding...")
                    warning.setForeground(QColor(self.secondary_colour_R, self.secondary_colour_G, self.secondary_colour_B, alpha=self.secondary_colour_A))
                    self.DialogBox.addItem(warning)
                else:
                    structure = ase.Atoms()
                    structure_dict = ast.literal_eval(self.structure[self.chosen_index])
                    for symbol, x, y, z in zip(structure_dict['symbol'], structure_dict['x_value'], structure_dict['y_value'], structure_dict['z_value']):
                        structure.append(ase.Atom(f"{symbol}", [float(x), float(y), float(z)]))
                    ase.visualize.view(structure)
            else:
                start_banner = QtWidgets.QListWidgetItem(f"\n!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING\n")
                start_banner.setForeground(QColor(self.secondary_colour_R, self.secondary_colour_G, self.secondary_colour_B, alpha=self.secondary_colour_A))
                self.DialogBox.addItem(start_banner)
                self.DialogBox.addItem("!!!WARNING!!! --> You are trying to visualize a benchmark complex. This is not supported --> Aborting Visualization")
                end_banner = QtWidgets.QListWidgetItem(f"\n!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!\n")
                end_banner.setForeground(QColor(self.secondary_colour_R, self.secondary_colour_G, self.secondary_colour_B, alpha=self.secondary_colour_A))
                self.DialogBox.addItem(end_banner)
                self.DialogBox.scrollToBottom()
        else:
            start_banner = QtWidgets.QListWidgetItem(f"\n!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING\n")
            start_banner.setForeground(QColor(self.secondary_colour_R, self.secondary_colour_G, self.secondary_colour_B, alpha=self.secondary_colour_A))
            self.DialogBox.addItem(start_banner)
            self.DialogBox.addItem("!!!WARNING!!! --> Please select a data point before attempting visualization --> Aborting Visualization")
            end_banner = QtWidgets.QListWidgetItem(f"\n!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!\n")
            end_banner.setForeground(QColor(self.secondary_colour_R, self.secondary_colour_G, self.secondary_colour_B, alpha=self.secondary_colour_A))
            self.DialogBox.addItem(end_banner)
            self.DialogBox.scrollToBottom()

    def Get_point(self, x_pos, y_pos):
        #####EXAMPLE#####
        smallest_distance = None
        dist_list = []
        # We loop through all the points to find which one is closest to our cursor
        for i, j, name in zip(self.x_axis_list, self.y_axis_list, self.Name):
            a = np.array((i, j))
            b = np.array((x_pos, y_pos))
            dist = np.linalg.norm(a - b)
            if smallest_distance is None:
                smallest_distance = dist
            if dist < smallest_distance:
                smallest_distance = dist
            else:
                pass
            dist_list.append(dist)

        index = dist_list.index(smallest_distance)
        self.chosen_index = index
        self.DialogBox.addItem(" ")
        start_banner = QtWidgets.QListWidgetItem(f"#####################---{self.Name[index]}---#####################")
        start_banner.setForeground(QColor(self.secondary_colour_R, self.secondary_colour_G, self.secondary_colour_B, alpha=self.secondary_colour_A))
        self.DialogBox.addItem(start_banner)
        self.DialogBox.addItem(" ")
        self.DialogBox.addItem(f"Complex Name: {self.Name[index]}")
        self.DialogBox.addItem(f"Crosshair X-value: {x_pos}")
        self.DialogBox.addItem(f"Crosshair Y-value: {y_pos}")
        self.DialogBox.addItem(f"Real X-value: {self.x_axis_list[index]}")
        self.DialogBox.addItem(f"Real Y-value: {self.y_axis_list[index]}")
        self.DialogBox.addItem(f"------------DATA------------\n")
        # self.DialogBox.addItem(f"------------DATA------------\n"f"{self.other_data[index]}")
        for line in self.other_data[index].split("\n"):
            self.DialogBox.addItem(str(line))
        self.DialogBox.addItem(f"------------END!------------\n")
        end_banner = QtWidgets.QListWidgetItem(f"##################################################################")
        end_banner.setForeground(QColor(self.secondary_colour_R, self.secondary_colour_G, self.secondary_colour_B, alpha=self.secondary_colour_A))
        self.DialogBox.addItem(end_banner)
        self.DialogBox.addItem(" ")
        self.DialogBox.scrollToBottom()

    def add_cross_hair(self):
        self.crosshair_v = pg.InfiniteLine(angle=90, movable=False, pen=pg.mkPen(color=self.secondary_colour_HEX, width=0.4))
        self.crosshair_h = pg.InfiniteLine(angle=0, movable=False, pen=pg.mkPen(color=self.secondary_colour_HEX, width=0.4))
        self.graphicsView.addItem(self.crosshair_v)
        self.graphicsView.addItem(self.crosshair_h)

    def clear(self):
        # todo: deep copy gives me issues here so I have resorted to making a new list
        plot_items = []
        for plot_item in self.graphicsView.plotItem.items:
            plot_items.append(plot_item)
        for plot_item in plot_items:
            self.graphicsView.removeItem(plot_item)
        self.add_cross_hair()

    def increment(self):
        self.clear()
        self.chosen_index = None
        if self.current_plot_number < len(self.plot_list) - 1:
            self.current_plot_number = self.current_plot_number + 1
        else:
            self.current_plot_number = 0
        self.plot_list[self.current_plot_number]()

    def decrement(self):
        self.clear()
        self.chosen_index = None
        if self.current_plot_number > 0:
            self.current_plot_number = self.current_plot_number - 1
        else:
            self.current_plot_number = len(self.plot_list) - 1
        self.plot_list[self.current_plot_number]()

    def mouseMovedEvent(self, event):
        pos = event[0]
        if self.graphicsView.sceneBoundingRect().contains(pos):
            mousePoint = self.graphicsView.getPlotItem().vb.mapSceneToView(pos)
            self.MOUSE_X = mousePoint.x()
            self.MOUSE_Y = mousePoint.y()

    def mouseClickedEvent(self):
        self.crosshair_v.setPos(self.MOUSE_X)
        self.crosshair_h.setPos(self.MOUSE_Y)
        self.Get_point(self.MOUSE_X, self.MOUSE_Y)

    def resizeEvent(self, event):
        self.resized.emit()
        return super(Window, self).resizeEvent(event)

    def someFunction(self):
        self.window_x = QWidget.size(self.centralWidget()).width()
        self.window_y = QWidget.size(self.centralWidget()).height()
        self.spacing = 20
        self.widget_ratio = 0.5
        self.graph_x_length = int((self.window_x - (1.5 * self.spacing)) * self.widget_ratio)
        self.button_height = int(3 * self.spacing)

        self.pushButton.setGeometry(QtCore.QRect(int(self.graph_x_length + (2 * self.spacing)), self.spacing, (self.window_x - ((3 * self.spacing) + self.graph_x_length)), self.button_height))

        self.pushButton_2a.setGeometry(
            QtCore.QRect(int(self.graph_x_length + (2 * self.spacing)), (2 * self.spacing) + self.button_height, int((self.window_x - ((5 * self.spacing) + self.graph_x_length)) / 3.0),
                         self.button_height))

        self.pushButton_2b.setGeometry(
            QtCore.QRect(int(self.graph_x_length + (2 * self.spacing)) + int((self.window_x - ((5 * self.spacing) + self.graph_x_length)) / 3.0) + self.spacing,
                         (2 * self.spacing) + self.button_height, int((self.window_x - ((5 * self.spacing) + self.graph_x_length)) / 3.0), self.button_height))

        self.pushButton_2c.setGeometry(
            QtCore.QRect(int(self.graph_x_length + (2 * self.spacing)) + 2 * int((self.window_x - ((5 * self.spacing) + self.graph_x_length)) / 3.0) + (2 * self.spacing),
                         (2 * self.spacing) + self.button_height, int((self.window_x - ((5 * self.spacing) + self.graph_x_length)) / 3.0), self.button_height))

        self.pushButton_3.setGeometry(
            QtCore.QRect(int(self.graph_x_length + (2 * self.spacing)), (3 * self.spacing) + (2 * self.button_height), (self.window_x - ((3 * self.spacing) + self.graph_x_length)),
                         self.button_height))

        self.DialogBox.setGeometry(
            QtCore.QRect(int(self.graph_x_length + (2 * self.spacing)), int((4 * self.spacing) + 3 * self.button_height), (self.window_x - ((3 * self.spacing) + self.graph_x_length)),
                         self.window_y - int((7 * self.spacing) + 3 * self.button_height)))

        self.graphicsView.setGeometry(QtCore.QRect(self.spacing, self.spacing, self.graph_x_length, self.window_y - (4 * self.spacing)))

    @staticmethod
    def hex_to_rgb(hex):
        new_hex = hex.strip("#")
        return tuple(int(new_hex[i:i + 2], 16) for i in (0, 2, 4))

    def view_chem_draw(self):
        # we need to make folder
        # todo if Timo ever gets a hold of this script I need to change os for pathlib
        assert self.current_plot_num is not None
        # Get the current working directory
        cwd = os.getcwd()
        if (self.downloadable_xyz_indexes is not None) and (self.downloadable_xyz_indexes != []):
            # Specify the name of the new directory
            new_directory = f"XYZ_FILES_PLOT_{self.current_plot_num}"
            # Create the new directory
            path = os.path.join(cwd, new_directory)
            os.mkdir(path)
            os.chdir(new_directory)

            print(f"Number of downloadable files = [{len(self.downloadable_xyz_indexes)}]")
            for index in self.downloadable_xyz_indexes:
                structure_dict = self.data[index]["Structure"]
                xyz_string = ""
                xyz_string = xyz_string + f"{len(structure_dict['symbol'])}" + "\n\n"
                for symbol, x, y, z in zip(structure_dict['symbol'], structure_dict['x_value'], structure_dict['y_value'], structure_dict['z_value']):
                    xyz_string = xyz_string + f"{symbol}" + "   " + str(float(x)) + "   " + str(float(y)) + "   " + str(float(z)) + "\n"
                with open(f"{index}" + ".xyz", 'w') as file:
                    file.write(xyz_string)
            start_banner = QtWidgets.QListWidgetItem(f"######################################################################")
            start_banner.setForeground(QColor(self.secondary_colour_R, self.secondary_colour_G, self.secondary_colour_B, alpha=self.secondary_colour_A))
            self.DialogBox.addItem(start_banner)
            self.DialogBox.addItem("!!!SUCCESS!!! -> XYZ files successfully downloaded")
            self.DialogBox.addItem(f"You can find them at:")
            self.DialogBox.addItem(f"{path}")
            end_banner = QtWidgets.QListWidgetItem(f"######################################################################")
            end_banner.setForeground(QColor(self.secondary_colour_R, self.secondary_colour_G, self.secondary_colour_B, alpha=self.secondary_colour_A))
            self.DialogBox.addItem(end_banner)
            self.DialogBox.scrollToBottom()
        else:
            start_banner = QtWidgets.QListWidgetItem(f"######################################################################")
            start_banner.setForeground(QColor(self.secondary_colour_R, self.secondary_colour_G, self.secondary_colour_B, alpha=self.secondary_colour_A))
            self.DialogBox.addItem(start_banner)
            self.DialogBox.addItem("!!!WARNING!!! -> XYZ downloading has not been enabled for this plot -> speak to Cian if you really need this :)")
            end_banner = QtWidgets.QListWidgetItem(f"######################################################################")
            end_banner.setForeground(QColor(self.secondary_colour_R, self.secondary_colour_G, self.secondary_colour_B, alpha=self.secondary_colour_A))
            self.DialogBox.addItem(end_banner)
            self.DialogBox.scrollToBottom()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    w = Window(parent=None, input_data=NBO_dict)
    w.show()
    sys.exit(app.exec_())
