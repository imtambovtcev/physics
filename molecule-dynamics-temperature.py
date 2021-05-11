import time
import moderngl as mgl
import numpy as np
from PyQt5 import QtOpenGL, QtWidgets, QtCore
import OpenGL.GL as gl
import sys
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import random

class Model():
    def __init__(self,size=30):
        self.size=int(size)
        self.set_on_lattice()
        self.t=1
        self.time=0.

    def set_on_lattice(self):
        nx=ny=int(np.ceil(np.sqrt(self.size)))
        rx=np.linspace(0,1-1/nx,nx)+0.5/nx
        print(f'{nx = }')
        print(f'{rx = }')
        ry = np.linspace(0, 1-1/ny, ny) + 0.5 / ny
        print(f'{ny = }')
        print(f'{ry = }')
        self.r=np.array(np.meshgrid(rx,ry)).reshape(2,-1).T
        self.r=self.r[:self.size,:]
        print(f'{self.r.shape = }')
        self.v = np.random.rand(self.size, 2)
        print(f'{self.r.shape = }')

    def set_rand(self):
        self.r=np.random.rand(self.size,2)
        self.v = np.random.rand(self.size, 2)
        print(f'{self.r.shape = }')

class Worker(QtCore.QRunnable):
    '''
    Worker thread
    '''

    @QtCore.pyqtSlot()
    def run(self):
        '''
        Your code goes in this function
        '''
        print("Thread start")
        time.sleep(5)
        print("Thread complete")

class RenderWidget(QtWidgets.QWidget):
    def __init__(self, *vargs, model=None, **kwargs):
        super().__init__(*vargs, **kwargs)
        self.record=False
        self.model = model
        self.setFixedSize(800, 800)

        # a figure instance to plot on
        self.figure = Figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `plot` method
        self.button = QtWidgets.QPushButton('Запись')
        self.button.clicked.connect(self.onrecord)

        # set the layout
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.button)
        self.setLayout(layout)
        self.updatePlot()

    def onrecord(self):
        if self.record:
            self.button.setText('Запись')
            self.record=False
        else:
            self.button.setText('Прервать запись')
            self.record = True
        print(self.figure)

    def updatePlot(self):
        ax = self.figure.add_subplot(111)
        ax.clear()
        ax.plot(*self.model.r.T, '*')
        ax.set_xlim([0,1])
        ax.set_ylim([0, 1])
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_title('t = {}'.format(self.model.time))
        self.canvas.draw()

######################################################################################################

class ParamWidget(QtWidgets.QWidget):
    def __init__(self, *vargs, model=None, **kwargs):
        super().__init__(*vargs, **kwargs)
        self.model = model

        self.randomWidget = QtWidgets.QPushButton("Случайное состояние")
        self.randomWidget.clicked.connect(self.onRandom)

        self.fmWidget = QtWidgets.QPushButton("Однородное состояние")
        self.fmWidget.clicked.connect(self.onFM)

        self.resetWidget = QtWidgets.QPushButton("Перезапустить")
        self.resetWidget.clicked.connect(self.onReset)

        self.tWidget = QtWidgets.QSlider(orientation=QtCore.Qt.Horizontal)
        self.tWidget.setMinimum(0)
        self.tWidget.setMaximum(500)
        self.tWidget.valueChanged.connect(self.ontValue)
        self.tWidget.setFixedWidth(200)


        vbox = QtWidgets.QVBoxLayout()
        # vbox.addWidget(QtWidgets.QLabel("Scene"))
        vbox.addWidget(self.resetWidget)
        vbox.addWidget(self.fmWidget)
        vbox.addWidget(self.randomWidget)
        vbox.addWidget(QtWidgets.QLabel("Температура"))
        vbox.addWidget(self.tWidget)

        vbox.addStretch(1)
        self.setLayout(vbox)

        self.reset_parameters()

        self.show()

    def reset_parameters(self):
        self.t=0

    def onReset(self, value):
        self.reset_parameters()

    def onFM(self, value):
        self.model.set_on_lattice()
        self.parent().updateEvent()

    def onRandom(self, value):
        self.model.set_rand()
        self.parent().updateEvent()

    def ontValue(self, value):
        self.model.t = value

    @property
    def t(self):
        return self.model.t

    @t.setter
    def t(self, value):
        self.model.t = value
        self.tWidget.setValue(int(self.t))


class MainWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.model=Model()
        # self.setFixedSize(self.WINDOW_SIZE[0], self.WINDOW_SIZE[1])
        self.move(QtWidgets.QDesktopWidget().rect().center() - self.rect().center())
        self.setWindowTitle("Моделирование газа")
        self.setWindowFlags(QtCore.Qt.Dialog)

        self.paramWidget = ParamWidget(model=self.model)
        self.renderWidget = RenderWidget(model=self.model)

        hbox = QtWidgets.QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(self.paramWidget)
        hbox.addWidget(self.renderWidget)
        self.setLayout(hbox)

        # self.paramWidget.onSceneChanged = self.renderWidget.load_scene

        self.show()

    def updateEvent(self):
        self.renderWidget.updatePlot()

    def keyPressEvent(self, event):
        # Quit when ESC is pressed
        if event.key() == QtCore.Qt.Key_Escape:
            QtCore.QCoreApplication.instance().quit()

    def keyReleaseEvent(self, event):
        pass

    def mouseMoveEvent(self, event):
        # self.wnd.mouse = (event.x(), event.y())
        pass

    def wheelEvent(self, event):
        # self.wnd.wheel += event.angleDelta().y()
        pass

##############################################################################################

if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    widget = MainWindow()
    sys.exit(app.exec_())