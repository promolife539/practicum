import sys
import re
import time
from controller import ParticleController, ParticlePlot

from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QIcon, QFont, QColor
from PyQt5.QtWidgets import (
    QApplication, QTextEdit, QLabel, QPushButton, QMessageBox, QDesktopWidget,
    QMainWindow, QSlider, QToolTip, QSizePolicy, QComboBox
)

class ParticleGUI(QMainWindow, ParticleController):
    handlers = { Qt.Key_Escape: lambda self: self.close() }
    qlwidth = 20
    qlheight = 24
    qtewidth = 55
    qteheight = 24

    def __init__(self):
        super().__init__(),
        QToolTip.setFont(QFont('Arial', 10))
        self.__init_main_window()
        self.__init_labels(left = 30, top = 400)
        self.__init_textboxes(left = 50, top = 400)
        self.__init_slider_particle_mass()
        self.__init_combobox_verlet()
        self.__init_button_add_particle()
        self.__init_button_stat()
        self.__init_plot()
        self.show()

    def __init_main_window(self):
        """ Init function for application's main window """
        self.resize(700, 400)
        self.setWindowTitle('Particles')
        frame = self.frameGeometry()
        frame.moveCenter(QDesktopWidget().availableGeometry().center())
        self.move(frame.topLeft())

    def __init_labels(self, left, top):
        s_lbl = QLabel("Velocity", self)
        s_lbl.setGeometry(left + __class__.qlwidth, top + __class__.qlheight,
            __class__.qtewidth, __class__.qlheight)
        u_lbl = QLabel(" u: ", self)
        u_lbl.setGeometry(left, top + 2 * __class__.qlheight, __class__.qlwidth,
            __class__.qlheight)
        v_lbl = QLabel(" v: ", self)
        v_lbl.setGeometry(left, top + 3.25 * __class__.qlheight, __class__.qlwidth,
            __class__.qlheight)
        p_lbl = QLabel("Position", self)
        p_lbl.setGeometry(left + 2 * __class__.qlwidth + __class__.qtewidth,
            top + __class__.qlheight, __class__.qtewidth, __class__.qlheight)
        x_lbl = QLabel(" x: ", self)
        x_lbl.setGeometry(left + __class__.qlwidth + __class__.qtewidth,
            top + 2 * __class__.qlheight, __class__.qlwidth, __class__.qlheight)
        y_lbl = QLabel(" y: ", self)
        y_lbl.setGeometry(left + __class__.qlwidth + __class__.qtewidth,
            top + 3.25 * __class__.qlheight, __class__.qlwidth, __class__.qlheight)
        e_lbl = QLabel("Optional", self)
        e_lbl.setGeometry(left + 3 * __class__.qlwidth + 2 * __class__.qtewidth,
            top + __class__.qlheight, __class__.qtewidth, __class__.qlheight)
        l_lbl = QLabel(" lifetime: ", self)
        l_lbl.setGeometry(left + 3 * __class__.qlwidth + 2 * __class__.qtewidth,
            top + 2 * __class__.qlheight, __class__.qtewidth, __class__.qlheight)
        m_lbl = QLabel(" mass: ", self)
        m_lbl.setGeometry(left + 3 * __class__.qlwidth + 2 * __class__.qtewidth,
            top + 3.25 * __class__.qlheight, __class__.qtewidth, __class__.qlheight)
        o_lbl = QLabel(" Methods: ", self)
        o_lbl.setGeometry(left + 5 * __class__.qlwidth + 4 * __class__.qtewidth,
            top + 3.25 * __class__.qlheight, 4 * __class__.qlwidth, __class__.qlheight)

    def __init_textboxes(self, left, top):
        self.u_tbx = QTextEdit(str(self.velocity["u"]), self)
        self.u_tbx.setGeometry(left, top + 2 * __class__.qlheight,
            __class__.qtewidth, __class__.qteheight)
        self.u_tbx.setTabChangesFocus(True)
        self.u_tbx.textChanged.connect(self.__u_tbx_changed)
        self.v_tbx = QTextEdit(str(self.velocity["v"]), self)
        self.v_tbx.setGeometry(left, top + 3.25 * __class__.qlheight,
            __class__.qtewidth, __class__.qteheight)
        self.v_tbx.setTabChangesFocus(True)
        self.v_tbx.textChanged.connect(self.__v_tbx_changed)
        self.x_tbx = QTextEdit(str(self.position["x"]), self)
        self.x_tbx.setGeometry(left + __class__.qlwidth + __class__.qtewidth,
            top + 2 * __class__.qteheight, __class__.qtewidth, __class__.qteheight)
        self.x_tbx.setTabChangesFocus(True)
        self.x_tbx.textChanged.connect(self.__x_tbx_changed)
        self.y_tbx = QTextEdit(str(self.position["y"]), self)
        self.y_tbx.setGeometry(left + __class__.qlwidth + __class__.qtewidth,
            top + 3.25 * __class__.qteheight, __class__.qtewidth, __class__.qteheight)
        self.y_tbx.setTabChangesFocus(True)
        self.y_tbx.textChanged.connect(self.__y_tbx_changed)
        self.l_tbx = QTextEdit(str(self.lifetime), self)
        self.l_tbx.setGeometry(left + 5 * __class__.qlwidth + 2 * __class__.qtewidth,
            top + 2 * __class__.qteheight, __class__.qtewidth, __class__.qteheight)
        self.l_tbx.setTabChangesFocus(True)
        self.l_tbx.textChanged.connect(self.__l_tbx_changed)

    def __init_slider_particle_mass(self):
        """Init function for slider that changes mass of the particle that will
        be created next"""
        sld = QSlider(Qt.Horizontal, self)
        sld.setFocusPolicy(Qt.NoFocus)
        sld.setGeometry(245, 475, 100, 30)
        sld.setMinimum(1)
        sld.setMaximum(1000)
        sld.setValue(self.mass / 5000)
        sld.valueChanged[int].connect(self.__sld_changed)

    def __init_combobox_verlet(self):
        """ Init combobox for selecting different Verlet implementation """
        self.cmb = QComboBox(self);
        self.cmb.setObjectName("cmb")
        self.cmb.setGeometry(420, 475, 250, 30)
        self.cmb.addItems(self.methods)
        self.cmb.currentIndexChanged.connect(self.__cmb_changed)

    def __init_button_add_particle(self):
        """ Init function for button that adds one particle to the plot """
        self.btn = QPushButton('Add particle', self)
        self.btn.setToolTip('Press this button to <b>Add a new particle</b>')
        self.btn.setGeometry(360, 430, 150, 40)
        self.btn.clicked.connect(self._ParticleController__add_particle)
        self.btn.setDisabled(False)

    def __init_plot(self):
        """ Init function for matplotlib plot """
        self.plot = ParticlePlot(self, 13.5, 8, 50, QSizePolicy.Fixed)
        self.plot.move(0, 0)
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.__draw_plot)
        self.timer.start(42)

    def __init_button_stat(self):
        """ Init function for button that shows benchmarks stats """
        self.sbtn = QPushButton('Yo, stats', self)
        self.sbtn.setToolTip('Press this button to <b>Show stats</b>')
        self.sbtn.setGeometry(520, 430, 150, 40)
        self.sbtn.clicked.connect(self.__sbtn_clicked)

    def __draw_plot(self):
        """ Function for particles rendering """
        self.plot.update_plot(self.particles, self.updaters[self.method])
        self.particles = list(filter(lambda p: p.time < p.death, self.particles))

    @staticmethod
    def __validate(textedit, pattern):
        if re.match(pattern, textedit.toPlainText()):
            textedit.setStyleSheet("QTextEdit {color: black}")
            return True
        else:
            textedit.setStyleSheet("QTextEdit {color: red}")
            return False

    def __u_tbx_changed(self):
        if __class__.__validate(self.u_tbx, r"^-?\d+(?:\.\d+)?$"):
            self.velocity["u"] = float(self.u_tbx.toPlainText())

    def __v_tbx_changed(self):
        if __class__.__validate(self.v_tbx, r"^-?\d+(?:\.\d+)?$"):
            self.velocity["v"] = float(self.v_tbx.toPlainText())

    def __x_tbx_changed(self):
        if __class__.__validate(self.x_tbx, r"^-?\d+(?:\.\d+)?$"):
            self.position["x"] = float(self.x_tbx.toPlainText())

    def __y_tbx_changed(self):
        if __class__.__validate(self.y_tbx, r"^-?\d+(?:\.\d+)?$"):
            self.position["y"] = float(self.y_tbx.toPlainText())

    def __l_tbx_changed(self):
        if __class__.__validate(self.l_tbx, r"^\d+?$"):
            self.lifetime = float(self.l_tbx.toPlainText())

    def __sld_changed(self, value):
        self.mass = value * 50000

    def __cmb_changed(self, value):
        self.method = value
        print("Current method: ", self.cmb.currentText())

    def __sbtn_clicked(self):
        print("============== Stats ============")
        print("Number of particles: {}".format(len(self.particles)))
        for i in range(3):
            start = time.time()
            self.updaters[i](self.particles)
            end = time.time()
            print("{} time = {}".format(self.updaters[i].__name__, end - start))
        print("============== Stats ============")

    def keyPressEvent(self, e):
        if e.key() in __class__.handlers.keys():
            __class__.handlers[e.key()](self)

    def closeEvent(self, event):
        reply = QMessageBox.question(self, 'Message', "Are you sure to quit?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        event.accept() if reply == QMessageBox.Yes else event.ignore()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    gui = ParticleGUI()
    sys.exit(app.exec())
