import os
import time
import moderngl as mgl
import numpy as np
from PyQt5 import QtOpenGL, QtWidgets, QtCore
import OpenGL.GL as gl
import sys
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib as mpl
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt

class Combo(QtWidgets.QComboBox):
    def __init__(self, dict):
        super(Combo, self).__init__()
        self.dict=dict
        self.keys=list(self.dict.keys())
        self.list=list(dict.values())
        print(f'{self.dict = }\n{self.keys = }\n{self.list = }\n')
        self.idx=0
        self.addItems(self.list)
        self.setCurrentIndex(self.idx)

    def currentElement(self):
        return self.keys[self.currentIndex()]

    def setCurrentElement(self,value):
        self.setCurrentIndex(self.keys.index(value))

    # @currentElement.setter
    # def currentEelment(self, value):
    #     self.model.t = value
    #     self.tWidget.setValue(int(self.t))




class State():
    def __init__(self):
        self.record=False
        self.pause=True
        self.record_n = 0
        self.update_request=False
        self.update_plot_request = False

class Model():
    def __init__(self,size=100):
        self._size = size
        self.t=250.
        self.time=0.
        self.data_save_rate=0.1
        self.last_data_update=0.
        self.interaction= None#'6_12'
        self.interaction_coef=-0.001
        self.boundary_list={'free':'Периодические','box':'Ящик'}
        self.boundary='free'
        self.r=None
        self.v=None
        self.contact=None
        self.vabs=None
        self.kinetic_energy=None
        self.set_on_lattice()
        self.updateInfo()
        self.data_names={'time':'$t$','mean V':'$< V >$','mean V squared':'$<V^2>$','mean r':'$<\Delta r>$'}
        self.data_captions = {'time': 'Время', 'mean V': 'Средняя скорость', 'mean V squared': 'Средний квадрат скорости','mean r': 'Среднее расстояние'}
        self.data_values = {'time': np.array([]), 'mean V': np.array([]),
                          'mean V squared': np.array([]),'mean r': np.array([])}
        self.data='mean V'

    @property
    def size(self):
        return self._size

    @size.setter
    def size(self,value):
        value=int(value)
        if value > self._size:
            self._size=value
        elif value<self._size:
            self._size = value

    def set_on_lattice(self):
        nx=ny=int(np.ceil(np.sqrt(self.size)))
        rx=np.linspace(0,1-1/nx,nx)+0.5/nx
        # print(f'{nx = }')
        # print(f'{rx = }')
        ry = np.linspace(0, 1-1/ny, ny) + 0.5 / ny
        # print(f'{ny = }')
        # print(f'{ry = }')
        self.r=np.array(np.meshgrid(rx,ry)).reshape(2,-1).T
        self.r=self.r[:self.size,:]
        self.v = 1.5*self.t/500*np.ones(self.r.shape)
        self.v[::2,0]*=-0.5
        #self.v = 2*(np.random.rand(self.size, 2)-0.5)
        self.updateInfo()

    def set_rand(self):
        self.r=np.random.rand(self.size,2)
        self.v = 2*(np.random.rand(self.size, 2)-0.5)
        self.updateInfo()

    def boundaryCheck(self):
        self.contact = self.r // 1
        if self.boundary == 'box':
            self.v-=2*self.v*np.abs(self.contact)
        else:
            self.r-=self.contact

    def run(self, dt):
        self.time+=dt
        if self.interaction == '6_12':
            all_r = np.concatenate([self.r, self.r + np.array([0, 1]), self.r + np.array([1, 0]),
                                    self.r - np.array([0, 1]), self.r - np.array([1, 0]),
                                    self.r + np.array([1, 1]), self.r + np.array([1, -1]),
                                    self.r - np.array([1, 1]), self.r - np.array([1, -1])])
            print(f'{all_r.shape = }')
            dr = np.array([np.subtract.outer(all_r[:, 0], all_r[:, 0]), np.subtract.outer(all_r[:, 1], all_r[:, 1])])
            dr = dr[:, :self.size, :]
            print(f'{dr.shape = }')
            print(f'{dr[:,1,:] = }')
            self.d = np.linalg.norm(dr, axis=0)
            di = np.divide(1., self.d)
            #d6 = np.power(di, 6)
            #d12 = np.power(d6, 2)
            #a = (d6 - d12)
            d3= np.power(di, 3)
            a=d3
            np.fill_diagonal(a, 0.)
            a1=np.sum(np.multiply(a,dr[0]),axis=1)
            a2 =np.sum(np.multiply(a,dr[0]),axis=1)
            print(f'{a.shape = } {a1.shape = } {dr.shape = }')
            a=np.array([a1,a2])
            a=self.interaction_coef*a.T
        else:
            a=0.

        self.v += dt*a
        self.r += dt * self.v

        self.boundaryCheck()


    def updateInfo(self):
        self.vabs=np.linalg.norm(self.v,axis=1)
        k=np.linalg.norm(self.vabs)
        self.kinetic_energy=k*k
        if self.interaction != '6_12':
            all_r = np.concatenate([self.r, self.r + np.array([0, 1]), self.r + np.array([1, 0]),
                                    self.r - np.array([0, 1]), self.r - np.array([1, 0]),
                                    self.r + np.array([1, 1]), self.r + np.array([1, -1]),
                                    self.r - np.array([1, 1]), self.r - np.array([1, -1])])
            dr = np.array([np.subtract.outer(all_r[:, 0], all_r[:, 0]), np.subtract.outer(all_r[:, 1], all_r[:, 1])])
            dr = dr[:, :self.size, :]
            self.d = np.linalg.norm(dr, axis=0)


    def updateData(self):
        if self.time-self.last_data_update>=self.data_save_rate:
            self.data_values['time']=np.append(self.data_values['time'],self.time)
            self.data_values['mean V']=np.append(self.data_values['mean V'],np.mean(self.vabs))
            self.data_values['mean V squared'] = np.append(self.data_values['mean V squared'], np.mean(np.multiply(self.vabs,self.vabs)))
            self.data_values['mean r'] = np.append(self.data_values['mean r'],
                                                           np.mean(self.d))
            self.last_data_update=self.time

class RenderWidget(QtWidgets.QWidget):
    def __init__(self, *vargs, model=None,state=None, **kwargs):
        super().__init__(*vargs, **kwargs)
        self.state=state
        self.model = model
        self.setFixedSize(800, 800)
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.add_subplot(111)
        self.norm = mpl.colors.Normalize(vmin=self.model.vabs.min(), vmax=self.model.vabs.max())
        self.cmap = mpl.cm.cool
        self.sm = plt.cm.ScalarMappable(cmap=self.cmap, norm=self.norm)
        #self.cbar=mpl.colorbar.ColorbarBase(self.ax, cmap=self.cmap, norm=self.norm)
        self.cbar=self.ax.figure.colorbar(self.sm)


        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)

        # Just some button connected to `plot` method
        self.button = QtWidgets.QPushButton('Запись')
        self.button.clicked.connect(self.onrecord)

        # set the layout
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.button)
        self.setLayout(layout)
        self.isfree=True
        self.updatePlot()
        self.timer = QtCore.QTimer()
        self.timer.start(int(1000/ 25))
        #self.timer.setInterval(1000)  # .5 seconds
        self.timer.timeout.connect(self.timePlotUpdate)

    def onrecord(self):
        if self.state.record:
            self.button.setText('Запись')
            self.state.record=False
        else:
            self.button.setText('Прервать запись')
            self.state.record = True
            self.save_directory = './record/'
            os.mkdir(self.save_directory)
        print(self.figure)

    def updatePlot(self,file=None):
        #print(f'plot update')
        if self.isfree:
            self.isfree=False
            self.ax.clear()
            self.norm = mpl.colors.Normalize(vmin=self.model.vabs.min(), vmax=self.model.vabs.max())
            self.ax.scatter(*self.model.r.T,c=self.model.vabs, norm=self.norm,cmap=self.cmap)
            self.sm.set_clim(self.model.vabs.min(), self.model.vabs.max())
            self.ax.set_xlim([0, 1])
            self.ax.set_ylim([0, 1])
            self.ax.set_title('t = {:10.4f}'.format(self.model.time))
            self.canvas.draw()
            self.isfree = True
            self.state.update_request=False
            if file is not None:
                self.state.record_n += 1
                self.figure.savefig(file)
            #print('Plot updated!')

    def timePlotUpdate(self):
        if not self.state.pause or self.state.update_plot_request:
            file=None if not self.state.record else self.save_directory+'{}'.format(self.state.record_n)
            self.updatePlot(file)



######################################################################################################

class ParamWidget(QtWidgets.QWidget):
    def __init__(self, *vargs, model=None,state=None, **kwargs):
        super().__init__(*vargs, **kwargs)
        self.model = model
        self.state=state

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
        self.pauseButton = QtWidgets.QPushButton('Запуск')
        self.pauseButton.clicked.connect(self.onPause)

        self.bWidget=Combo(self.model.boundary_list)
        self.bWidget.currentIndexChanged.connect(self.onbValue)

        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.resetWidget)
        vbox.addWidget(self.fmWidget)
        vbox.addWidget(self.randomWidget)
        vbox.addWidget(QtWidgets.QLabel("Температура"))
        vbox.addWidget(self.tWidget)
        vbox.addWidget(self.bWidget)
        vbox.addWidget(self.pauseButton)
        vbox.addStretch(1)
        self.setLayout(vbox)
        self.reset_parameters()
        self.show()

        self.isfree = True
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.updateState)

    def updateState(self):
        if self.isfree:
            self.isfree = False
            for _ in range(10):
                self.model.run(0.001)
            self.model.updateInfo()
            self.model.updateData()
            self.isfree = True
            #print('State updated!')

    def reset_parameters(self):
        self.model.time=0.
        self.t=250

    def onReset(self, value):
        self.reset_parameters()

    def onFM(self, value):
        self.model.set_on_lattice()
        self.state.update_request=True

    def onRandom(self, value):
        self.model.set_rand()
        self.state.update_request = True

    def ontValue(self, value):
        self.model.t = value

    def onbValue(self):
        self.model.boundary=self.bWidget.currentElement()
        print(f'{self.model.boundary = }')

    def onPause(self):
        if self.state.pause:
            self.timer.start(int(1000 / 10))
            self.state.pause = False
            self.pauseButton.setText('Пауза')
        else:
            self.timer.stop()
            self.state.pause = True
            self.pauseButton.setText('Запуск')

    @property
    def t(self):
        return self.model.t

    @t.setter
    def t(self, value):
        self.model.t = value
        self.tWidget.setValue(int(self.t))

    @property
    def boundary(self):
        return self.model.boundary

    @boundary.setter
    def boundary(self, value):
        self.model.boundary = value
        self.tWidget.setValue(value)

class GraphWindow(QtWidgets.QWidget):
    def __init__(self, parent,model=None,state=None):
        super().__init__(parent)
        self.model=model
        self.state=state
        self.current_plot='mean V'
        self.setWindowTitle("Моделирование газа")
        self.setWindowFlags(QtCore.Qt.Dialog)
        self.setFixedSize(800, 800)
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.add_subplot(111)

        self.combo=Combo(self.model.data_captions)
        self.combo.setCurrentElement(self.current_plot)
        self.combo.currentIndexChanged.connect(self.onchange)

        # set the layout
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.combo)
        self.setLayout(layout)
        self.isfree=True
        self.updatePlot()
        self.timer = QtCore.QTimer()
        self.timer.start(int(1000/ 25))
        self.timer.timeout.connect(self.timePlotUpdate)
        self.show()

    def onchange(self):
        self.current_plot=self.combo.currentElement()
        self.updatePlot()

    def updatePlot(self):
        if self.isfree:
            self.isfree=False
            self.ax.clear()
            self.ax.plot(self.model.data_values['time'],self.model.data_values[self.current_plot])
            self.ax.set_xlabel('$t$')
            self.ax.set_ylabel(self.model.data_names[self.current_plot])
            self.canvas.draw()
            self.isfree = True
            self.state.update_plot_request=False

    def timePlotUpdate(self):
        if not self.state.pause or self.state.update_request:
            self.updatePlot()


class MainWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.model=Model()
        self.state=State()
        # self.setFixedSize(self.WINDOW_SIZE[0], self.WINDOW_SIZE[1])
        self.move(QtWidgets.QDesktopWidget().rect().center() - self.rect().center())
        self.setWindowTitle("Моделирование газа")
        self.setWindowFlags(QtCore.Qt.Dialog)

        self.paramWidget = ParamWidget(model=self.model,state=self.state)
        self.renderWidget = RenderWidget(model=self.model,state=self.state)
        self.graphWindow = GraphWindow(self,model=self.model,state=self.state)

        hbox = QtWidgets.QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(self.paramWidget)
        hbox.addWidget(self.renderWidget)
        self.setLayout(hbox)
        # self.paramWidget.onSceneChanged = self.renderWidget.load_scene
        self.show()

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