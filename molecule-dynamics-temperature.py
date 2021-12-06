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
from scipy.stats import maxwell


# класс для выпадающей менюшки

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


# класс для состояния записи

class State():
    def __init__(self):
        self.record=False
        self.pause=True
        self.record_n = 0
        self.update_request=False
        self.update_plot_request = False

# основной класс для хранения данных

class Model():
    def __init__(self,size=100):
        self.t=250. # температура
        self.time=0. #время
        self.data_save_rate=0.01 # стандартная частота сохранения
        self.data_save_rate_coeff = 1. # во сколько раз частота сохранения отличается от исходной
        self.last_data_update=0. # когда было последнее обновление данных (в времени моделирования)
        self.reservoir_energy=0. # энергия резервуара (для моелирования взаимодействия с прилипанем)
        self.interaction_coef1 = 0.0001 # коэффициетны в кулоновском взаимодействии
        self.interaction_coef=0.001
        self.interaction_dict = {'No':'Нет','hard_balls': 'Твердые сферы', '-2':'Кулоновское взаимодействие'}  # описание всего, что отображается в меню взаимодействий
        self.interaction = 'No' # текущий выбранный пункт
        self.boundary_dict={'free':'Периодические','box':'Коробка','box_w':'Коробка с прилипанием'} # взаимодействия с поверхностью
        self.boundary='free'
        self.r=None # координаты частиц
        self.v=None # скорости частиц
        self.contact=None # массив для хранения частиц, ударившихся о стенку за последний тик
        self.vabs=None # модули скоростей
        self.kinetic_energy=None # кинетические энергии
        self.size = size # количество частиц
        self.updateInfo() #при инициализации сразу обновляется состояние
        self.data_names={'time':'$t$','mean V':'$< V >$','mean V squared':'$<V^2>$','mean r':'$<\Delta r>$','V dist':'$|V|$','Vx dist':'$V_x$','r dist':'$\Delta r$',
                         'first_moment': 'Первый момент','second_moment': 'Второй момент','third_moment': 'Третий момент','forth_moment': 'Четвертый момент','fifth_moment': 'Пятый момент'} # меню графиков
        self.data_captions = {'time': 'Время', 'mean V': 'Средняя скорость', 'mean V squared': 'Средний квадрат скорости',
                              'mean r': 'Среднее расстояние','V dist':'Распределение модулей скоростей','Vx dist':'Распределение проекций скоростей','r dist':'Распределение расстояний',
                              'first_moment': 'Первый момент','second_moment': 'Второй момент','third_moment': 'Третий момент','forth_moment': 'Четвертый момент','fifth_moment': 'Пятый момент'}
        self.data_values = {'time': np.array([]), 'mean V': np.array([]),
                          'mean V squared': np.array([]),'mean r': np.array([]),'V dist': np.array([]),'Vx dist': np.array([]),'r dist': np.array([]),
                            'first_moment':np.array([]),'second_moment': np.array([]),'third_moment': np.array([]),'forth_moment': np.array([]),'fifth_moment': np.array([])}
        self.data='mean V'

#размер удобнее устанавливать через setter getter
    @property
    def size(self):
        return self._size

    @size.setter
    def size(self,value):
        value=int(value)
        self._size = value
        self.set_on_lattice()
        print(f'{self.size = }')

# установка частиц на решетке
    def set_on_lattice(self):
        nx=ny=int(np.ceil(np.sqrt(self.size))) # задание размера решетки
        rx=np.linspace(0,1-1/nx,nx)+0.5/nx # задание позиций частиц на ней
        # print(f'{nx = }')
        # print(f'{rx = }')
        ry = np.linspace(0, 1-1/ny, ny) + 0.5 / ny
        # print(f'{ny = }')
        # print(f'{ry = }')
        self.r=np.array(np.meshgrid(rx,ry)).reshape(2,-1).T #создание частиц на решетке (их сейчас больше чем надо)
        self.r=self.r[:self.size,:] #ообрезание лишних
        self.v = 1.5*np.sqrt(self.t/500)*np.ones(self.r.shape) #зададим скорости
        self.v[::2,0]*=-0.2
        self.v[::2,1]*=0.5#
        self.v[:, 0] -= self.v[:, 0].mean()
        self.v[:, 1] -= self.v[:, 1].mean()
        #self.v = 2*(np.random.rand(self.size, 2)-0.5)
        self.updateInfo()

    def set_rand(self):
        self.r=np.random.rand(self.size,2)
        self.v =np.random.rand(2,self.size)
        print(f'{self.v[0].shape = }')
        self.v = np.array([np.multiply(self.v[0],np.cos(2*np.pi*self.v[1])),np.multiply(self.v[0],np.sin(2*np.pi*self.v[1]))]).T
        self.v *= 2*np.sqrt(self.t/85)
        self.v[:,0]-=self.v[:,0].mean()
        self.v[:, 1] -= self.v[:, 1].mean()
        print(f'{self.v.shape = }')
        self.updateInfo()

    def boundaryCheck(self):
        self.contact = self.r // 1
        if self.boundary == 'box':
            self.v-=2*self.v*np.abs(self.contact)
        elif self.boundary == 'box_w':
            is_out=(np.sum((self.contact!=0),axis=1))!=0
            is_out_n=np.sum(is_out)
            if is_out_n>0:
                #print(f'{np.mean(np.abs(self.v)) = }')
                #print(f'{is_out = }\n{is_out_n = }')
                #print(f'{np.linalg.norm(self.v[is_out,:],axis=1).shape = }')
                #print(f'{np.linalg.norm(self.v[is_out,:],axis=1) = }\n{np.linalg.norm(self.v[is_out,:],axis=1)*np.linalg.norm(self.v[is_out,:],axis=1) = }')
                energy_out=np.sum(np.linalg.norm(self.v[is_out,:],axis=1)*np.linalg.norm(self.v[is_out,:],axis=1))+self.reservoir_energy
                #print(f'{energy_out/(is_out_n+1) = }')
                #abs_v=1+np.random.standard_normal([is_out_n+1])
                abs_v=maxwell.rvs(size=is_out_n+1)
                abs_v*=np.sqrt(energy_out/np.sum(abs_v*abs_v))
                #print(f'{abs_v*abs_v = }')
                self.reservoir_energy=abs_v[-1]*abs_v[-1]
                abs_v=abs_v[:-1]
                #print(f'{abs_v = }')
                #print(f'{self.reservoir_energy = }')
                ang=2*np.pi*np.random.rand(is_out_n)
                new_velocities=np.array([np.multiply(abs_v,np.cos(ang)),np.multiply(abs_v,np.sin(ang))]).T
                #print(f'{new_velocities = }')
                #print(f'{self.contact[is_out] = }\n{self.contact[is_out] != 0 = }\n{self.contact[self.contact!=0] = }')
                new_velocities[self.contact[is_out] != 0] = -self.contact[self.contact!=0]*np.abs(new_velocities[self.contact[is_out] != 0])
                #print(f'{self.v[is_out,:].shape = }')
                self.v[is_out,:]=new_velocities
                #self.v[self.contact!=0]=-self.contact*self.v[self.contact!=0]
                #print(f'{new_velocities = }\n{energy_rand = }')
            #self.v -= 2 * self.v * np.abs(self.contact)
        else:
            self.r-=self.contact

    def run(self, dt):
        self.time+=dt
        if self.interaction == 'hard_balls':
            dm=distance_matrix(self.r,self.r)
            #print(f'{dm = }\n{dm.shape = }')
            interaction_distance=0.015
            #print(f'{dm.mean() = }')
            #print(f'{dm < interaction_distance = }')
            #print(f'{np.argwhere(dm < interaction_distance) = }')
            for i,j in np.argwhere(dm < interaction_distance):
                if i>j:
                    absv1=np.linalg.norm(self.v[i])
                    absv2=np.linalg.norm(self.v[j])
                    theta1=np.arctan2(self.v[i][1],self.v[i][0])
                    theta2 = np.arctan2(self.v[j][1], self.v[j][0])
                    #dr=self.r[i]-self.r[j]
                    phi =  theta2-theta1
                    #phi=np.arctan2(dr[1],dr[0])
                    #print(f'{theta1 = } {theta2 = } {phi = }')
                    self.v[i]=[absv2*np.cos(theta2-phi)*np.cos(phi)+absv1*np.sin(theta1-phi)*np.cos(phi+0.5*np.pi),
                               absv2*np.cos(theta2-phi)*np.sin(phi)+absv1*np.sin(theta1-phi)*np.sin(phi+0.5*np.pi)]
                    self.v[j]=[absv1 * np.cos(theta1 - phi) * np.cos(phi) + absv2 * np.sin(theta2 - phi) * np.cos(phi + 0.5 * np.pi),
                               absv1 * np.cos(theta1 - phi) * np.sin(phi) + absv2 * np.sin(theta2 - phi) * np.sin(
                                   phi + 0.5 * np.pi)]
                    #print(f'after: {self.v[i] = } {self.v[j] = }')

            a = 0.
        elif self.interaction == '-2':
            if self.boundary == 'free':
                all_r = np.concatenate([self.r, self.r + np.array([0, 1]), self.r + np.array([1, 0]),
                                        self.r - np.array([0, 1]), self.r - np.array([1, 0]),
                                        self.r + np.array([1, 1]), self.r + np.array([1, -1]),
                                        self.r - np.array([1, 1]), self.r - np.array([1, -1])])
                #print(f'{all_r.shape = }')
                dr = np.array([np.subtract.outer(all_r[:, 0], all_r[:, 0]), np.subtract.outer(all_r[:, 1], all_r[:, 1])])
                dr = dr[:, :self.size, :]
                #print(f'{dr.shape = }')
                #print(f'{dr[:,1,:] = }')
                self.d = np.linalg.norm(dr, axis=0)
                di = np.divide(1., self.d)
                #d6 = np.power(di, 6)
                #d12 = np.power(d6, 2)
                #a = (d6 - d12)
                d3= np.power(di, 3)
                a=d3
                np.fill_diagonal(a, 0.)
                a *= self.interaction_coef1
                a1=np.sum(np.multiply(a,dr[0]),axis=1)
                a2 =np.sum(np.multiply(a,dr[0]),axis=1)
                #print(f'{a.shape = } {a1.shape = } {dr.shape = }')
                a=np.array([a1,a2])
                a=self.interaction_coef*a.T
            else:
                dr = np.array(
                    [np.subtract.outer(self.r[:, 0], self.r[:, 0]), np.subtract.outer(self.r[:, 1], self.r[:, 1])])
                self.d = np.linalg.norm(dr, axis=0)
                di = np.divide(1., self.d)
                d3 = np.power(di, 3)
                a = d3
                np.fill_diagonal(a, 0.)
                a*=self.interaction_coef1
                #print(f'{a.max() = }')
                a1 = np.sum(np.multiply(a, dr[0]), axis=1)
                a2 = np.sum(np.multiply(a, dr[1]), axis=1)
                #print(f'{a.shape = } {a1.shape = } {dr.shape = }')
                a = np.array([a1, a2])
                a = self.interaction_coef * a.T
                #print(f'{a.max() = }')
        else:
            a=0.

        self.v += dt*a
        self.r += dt * self.v

        self.boundaryCheck()


    def updateInfo(self):
        self.vabs=np.linalg.norm(self.v,axis=1)
        k=np.linalg.norm(self.vabs)
        self.kinetic_energy=k*k
        if not self.interaction == '-2':
            if self.boundary == 'free':
                all_r = np.concatenate([self.r, self.r + np.array([0, 1]), self.r + np.array([1, 0]),
                                        self.r - np.array([0, 1]), self.r - np.array([1, 0]),
                                        self.r + np.array([1, 1]), self.r + np.array([1, -1]),
                                        self.r - np.array([1, 1]), self.r - np.array([1, -1])])
                dr = np.array([np.subtract.outer(all_r[:, 0], all_r[:, 0]), np.subtract.outer(all_r[:, 1], all_r[:, 1])])
                dr = dr[:, :self.size, :]
                self.d = np.linalg.norm(dr, axis=0)
            else:
                dr = np.array(
                    [np.subtract.outer(self.r[:, 0], self.r[:, 0]), np.subtract.outer(self.r[:, 1], self.r[:, 1])])
                self.d = np.linalg.norm(dr, axis=0)


    def updateData(self):
        if self.time-self.last_data_update>=self.data_save_rate*self.data_save_rate_coeff:
            self.data_values['time']=np.append(self.data_values['time'],self.time)
            self.data_values['mean V']=np.append(self.data_values['mean V'],np.mean(self.vabs))
            self.data_values['mean V squared'] = np.append(self.data_values['mean V squared'], np.mean(np.multiply(self.vabs,self.vabs)))
            self.data_values['mean r'] = np.append(self.data_values['mean r'], np.mean(self.d))
            self.data_values['first_moment'] = np.append(self.data_values['first_moment'], np.mean(self.v[:,0]))
            self.data_values['second_moment'] = np.append(self.data_values['second_moment'], np.mean(self.v[:, 0]**2))
            self.data_values['third_moment'] = np.append(self.data_values['third_moment'], np.mean(self.v[:, 0]**3))
            self.data_values['forth_moment'] = np.append(self.data_values['forth_moment'], np.mean(self.v[:, 0]**4))
            self.data_values['fifth_moment'] = np.append(self.data_values['fifth_moment'], np.mean(self.v[:, 0]**5))
            self.data_values['V dist'] = np.copy(self.vabs)
            self.data_values['Vx dist'] = self.v[:,0]
            self.data_values['r dist'] = np.copy(self.d).reshape(-1)
            self.last_data_update=self.time

            if len(self.data_values['time'])>200:
                for key in self.data_values.keys():
                    self.data_values[key] = self.data_values[key][::2]
                self.data_save_rate_coeff*=2
                #print(f'{self.data_save_rate_coeff = }\n{self.data_save_rate*self.data_save_rate_coeff = }')

    def clearData(self):
        self.data_values = {'time': np.array([]), 'mean V': np.array([]),
                            'mean V squared': np.array([]), 'mean r': np.array([]), 'V dist': np.array([]), 'Vx dist': np.array([]),'r dist': np.array([]),
                            'first_moment': np.array([]),'second_moment': np.array([]),'third_moment': np.array([]),'forth_moment': np.array([]),'fifth_moment': np.array([])}
        self.data_save_rate_coeff=1.
        self.last_data_update=0.

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
        self.cbar.ax.set_title('Скорость', fontsize=10)


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
            try:
                os.mkdir(self.save_directory)
            except:
                pass
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
        if (not self.state.pause) or self.state.update_request:
            file=None if not self.state.record else self.save_directory+'{:10.4f}.png'.format(self.model.time)
            self.updatePlot(file)



class ParamWidget(QtWidgets.QWidget):
    def __init__(self, *vargs, model=None,state=None, **kwargs):
        super().__init__(*vargs, **kwargs)
        self.model = model
        self.state=state

        self.isfree = False
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.updateState)

        self.randomWidget = QtWidgets.QPushButton("Случайное состояние")
        self.randomWidget.clicked.connect(self.onRandom)

        self.fmWidget = QtWidgets.QPushButton("Однородное состояние")
        self.fmWidget.clicked.connect(self.onFM)

        self.resetWidget = QtWidgets.QPushButton("Перезапустить")
        self.resetWidget.clicked.connect(self.onReset)

        self.nWidget = QtWidgets.QSlider(orientation=QtCore.Qt.Horizontal)
        self.nWidget.setMinimum(4)
        self.nWidget.setMaximum(500)
        self.nWidget.valueChanged.connect(self.onnValue)
        self.nWidget.setFixedWidth(200)
        self.tWidget = QtWidgets.QSlider(orientation=QtCore.Qt.Horizontal)
        self.tWidget.setMinimum(0)
        self.tWidget.setMaximum(500)
        self.tWidget.valueChanged.connect(self.ontValue)
        self.tWidget.setFixedWidth(200)
        self.sWidget = QtWidgets.QSlider(orientation=QtCore.Qt.Horizontal)
        self.sWidget.setMinimum(0)
        self.sWidget.setMaximum(20)
        self.sWidget.valueChanged.connect(self.onsValue)
        self.sWidget.setFixedWidth(200)

        self.pauseButton = QtWidgets.QPushButton('Запуск')
        self.pauseButton.clicked.connect(self.onPause)

        self.bWidget=Combo(self.model.boundary_dict)
        self.bWidget.currentIndexChanged.connect(self.onbValue)
        self.iWidget = Combo(self.model.interaction_dict)
        self.iWidget.currentIndexChanged.connect(self.oniValue)

        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.resetWidget)
        vbox.addWidget(self.fmWidget)
        vbox.addWidget(self.randomWidget)
        vbox.addWidget(QtWidgets.QLabel("Количество частиц"))
        vbox.addWidget(self.nWidget)
        vbox.addWidget(QtWidgets.QLabel("Температура"))
        vbox.addWidget(self.tWidget)
        vbox.addWidget(QtWidgets.QLabel("Граничные условия"))
        vbox.addWidget(self.bWidget)
        vbox.addWidget(QtWidgets.QLabel("Взаимодействие"))
        vbox.addWidget(self.iWidget)
        vbox.addWidget(QtWidgets.QLabel("Скорость воспроизведения"))
        vbox.addWidget(self.sWidget)
        vbox.addWidget(self.pauseButton)
        vbox.addStretch(1)
        self.setLayout(vbox)
        self.reset_parameters()
        self.show()
        self.isfree = True

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
        self.n=self.model.size
        self.t=250
        self.s = 0
        self.pause=True
        self.state.update_request = True
        self.state.update_plots_request = True

    def onReset(self, value):
        self.reset_parameters()
        self.model.clearData()

    def onFM(self, value):
        self.model.set_on_lattice()
        self.state.update_request=True
        self.state.update_plots_request = True

    def onRandom(self, value):
        self.model.set_rand()
        self.state.update_request = True
        self.state.update_plots_request = True

    def onnValue(self, value):
        self.n = value
        self.state.update_request = True
        self.state.update_plots_request = True

    def ontValue(self, value):
        self.t = value

    def onsValue(self, value):
        self.s = value

    def onbValue(self):
        self.model.boundary=self.bWidget.currentElement()
        print(f'{self.model.boundary = }')

    def oniValue(self):
        self.model.interaction=self.iWidget.currentElement()
        print(f'{self.model.interaction = }')

    def onPause(self):
        self.pause=not(self.state.pause)

    @property
    def pause(self):
        return self.state.pause

    @pause.setter
    def pause(self,value):
        self.state.pause=value
        if self.state.pause:
            self.timer.stop()
            self.state.pause = True
            self.pauseButton.setText('Запуск')
        else:
            self.timer.start(int(self.s)*10)
            self.state.pause = False
            self.pauseButton.setText('Пауза')

    @property
    def n(self):
        return self.model.size

    @n.setter
    def n(self, value):
        self.model.size = value
        self.nWidget.setValue(int(self.n))
        print(f'{self.n = }')

    @property
    def t(self):
        return self.model.t

    @t.setter
    def t(self, value):
        self.model.t = value
        self.tWidget.setValue(int(self.t))
        print(f'{self.t = }')

    @property
    def s(self):
        return self._s

    @s.setter
    def s(self, value):
        self._s = value
        self.sWidget.setValue(int(self._s))
        self.timer.setInterval(int(self._s)*10)
        print(f'{self.s = }')

    @property
    def boundary(self):
        return self.model.boundary

    @boundary.setter
    def boundary(self, value):
        self.model.boundary = value
        self.tWidget.setValue(value)
        print(f'{self.model.boundary = }')

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
            if self.current_plot == 'V dist':
                self.ax.hist(self.model.data_values[self.current_plot],bins=20)
                self.ax.set_xlabel('$|V|$')
                self.ax.set_ylabel('Количество частиц')
            elif self.current_plot == 'Vx dist':
                self.ax.hist(self.model.data_values[self.current_plot],bins=20)
                self.ax.set_xlabel('$V_x$')
                self.ax.set_ylabel('Количество частиц')
            elif self.current_plot == 'r dist':
                self.ax.hist(self.model.data_values[self.current_plot],bins=20)
                self.ax.set_xlabel('$r$')
                self.ax.set_ylabel('Количество частиц')
            else:
                self.ax.plot(self.model.data_values['time'],self.model.data_values[self.current_plot])
                self.ax.set_xlabel('$t$')
                self.ax.set_ylabel(self.model.data_names[self.current_plot])
            self.canvas.draw()
            self.isfree = True
            self.state.update_plots_request=False

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