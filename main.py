import math
import sys
import numpy as np
from PyQt6.QtGui import QPalette, QColor
from PyQt6.QtWidgets import QMainWindow, QMessageBox, QVBoxLayout, QTableWidget, QTableWidgetItem, QApplication
from PyQt6.uic import loadUi
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
from matplotlib.patches import Polygon as MPolygon, Circle as MCircle


def point_A_prime_function(phi, fcc_prime, A_g, A_s, f_y):
    try:
        return phi * 0.8 * (0.85 * fcc_prime * (A_g - A_s) + f_y * A_s)
    except ZeroDivisionError:
        print('Can\'t divide by zero.')


def A_function(b, E_c, E_2, fc_prime, eccu, c):
    try:
        return -b * ((E_c - E_2) ** 2 / (12 * fc_prime)) * (eccu / c) ** 2
    except ZeroDivisionError:
        print('Can\'t divide by zero.')


def B_function(b, E_c, E_2, eccu, c):
    try:
        return b * (E_c - E_2) / 2 * (eccu / c)
    except ZeroDivisionError:
        print('Can\'t divide by zero.')


def C_function(b, fc_prime):
    return -b * fc_prime


def D_function(b, c, fc_prime, E_2, eccu):
    try:
        return b * c * (fc_prime + E_2 / 2 * eccu)
    except ZeroDivisionError:
        print('Can\'t divide by zero.')
        return None
    else:  # inserted
        pass


def E_function(b, E_c, E_2, fc_prime, eccu, c):
    try:
        return -b * (E_c - E_2) ** 2 / (16 * fc_prime) * (eccu / c) ** 2
    except ZeroDivisionError:
        print('Can\'t divide by zero.')


def F_function(b, h, c, E_c, E_2, fc_prime, eccu):
    try:
        return b * ((c - h / 2) * ((E_c - E_2) ** 2 / (12 * fc_prime)) * (eccu / c) ** 2 + (E_c - E_2) / 3 * (eccu / c))
    except ZeroDivisionError:
        print('Can\'t divide by zero.')


def G_function(b, h, c, E_c, E_2, fc_prime, eccu):
    try:
        return -(b / 2 * fc_prime + b * (c - h / 2) * ((E_c - E_2) / 2) * (eccu / c))
    except ZeroDivisionError:
        print('Can\'t divide by zero.')


def H_function(b, h, c, fc_prime):
    try:
        return b * fc_prime * (c - h / 2)
    except ZeroDivisionError:
        print('Can\'t divide by zero.')


def I_function(b, h, c, E_2, fc_prime, eccu):
    try:
        return b * c ** 2 * fc_prime / 2 - b * c * fc_prime * (c - h / 2) + b * c ** 2 * E_2 * eccu / 3 - b * c * E_2 / 2 * (c - h / 2) * eccu
    except ZeroDivisionError:
        print('Can\'t divide by zero.')


def Pn_function(A, B, C, D, A_s, f_s, yt):
    try:
        return A * yt ** 3 + B * yt ** 2 + C * yt + D + A_s * f_s
    except ZeroDivisionError:
        print('Can\'t divide by zero.')


def Mn_function(E, F, G, H, I, side_bars_forces, yt):
    try:
        return E * yt ** 4 + F * yt ** 3 + G * yt ** 2 + H * yt + I + side_bars_forces
    except ZeroDivisionError:
        print('Can\'t divide by zero.')


def steel_strain(c, cover, spacing, layer, eccu, E_s, f_y):
    # Calculate yield strain once
    yield_strain = f_y / E_s

    # Calculate actual strain based on the input parameters
    actual_strain = eccu * ((c - cover - spacing * layer) / c)

    # Return the minimum of actual strain and yield strain
    return min(abs(actual_strain), yield_strain) * (1 if actual_strain >= 0 else -1)


class MainApp(QMainWindow):
    def __init__(self, parent=None):
        super(MainApp, self).__init__(parent)
        loadUi('gui.ui', self)
        self.setWindowTitle('FRP-confined RC column.')
        palette = self.tableWidget.palette()
        palette.setColor(QPalette.ColorRole.AlternateBase, QColor('#1f1f1f'))
        self.tableWidget.setPalette(palette)
        palette_2 = self.tableWidget_2.palette()
        palette_2.setColor(QPalette.ColorRole.AlternateBase, QColor('#1f1f1f'))
        self.tableWidget_2.setPalette(palette_2)
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        self.setup_canvas()
        self.points = []
        self.exposure = 1.0
        self.phi = 1.0
        self.t_f = 1.0
        self.Ef = 1.0
        self.eps_fu_manu = 1.0
        self.f_fu_manu = 1.0
        self.f_fu = 1.0
        self.E_s = 1.0
        self.f_y = 1.0
        self.n_plies = 1.0
        self.psi_f = 1.0
        self.frp_mat = 1.0
        self.area_of_main_rebar = 1.0
        self.area_of_sec_rebar = 1.0
        self.area_of_side_rebar = 1.0
        self.fc_prime = 1.0
        self.number_of_main_rebar = 1
        self.number_of_side_rebar = 1
        self.number_of_sec_rebar = 1
        self.dia_of_side_rebar = 1.0
        self.dia_of_sec_rebar = 1.0
        self.dia_of_main_rebar = 1.0
        self.r_c = 1.0
        self.b = 1.0
        self.h = 1.0
        self.A_g = 1.0
        self.CE = 1.0
        self.A_s = 1.0
        self.epsilon_fu = 1.0
        self.epsilon_fe = 1.0
        self.diagonal = 1.0
        self.f_l = 1.0
        self.A_c = 1.0
        self.rho_g = 1.0
        self.A_e = 1.0
        self.k_a = 1.0
        self.k_b = 1.0
        self.epsilon_c_prime = 1.0
        self.epsilon_f_y = 1.0
        self.k_e = 1.0
        self.d = 1.0
        self.c_point_C = 1.0
        self.E_c = 1.0
        self.diameter = 1.0
        self.epsilon_fe_POINT_A_prime = 1.0
        self.epsilon_fe_POINT_B = 1.0
        self.epsilon_fe_POINT_C = 1.0
        self.epsilon_fe_POINT_D = 1.0
        self.f_l_POINT_A_prime = 1.0
        self.f_l_POINT_B = 1.0
        self.f_l_POINT_C = 1.0
        self.f_l_POINT_D = 1.0
        self.epsilon_ccu_POINT_A_prime = 1.0
        self.epsilon_ccu_POINT_B = 1.0
        self.epsilon_ccu_POINT_C = 1.0
        self.epsilon_ccu_POINT_D = 1.0
        self.bars_spacing = 1.0
        self.cover = 1.0
        self.E_2_A = 1.0
        self.E_2_B = 1.0
        self.E_2_C = 1.0
        self.E_2_D = 1.0
        self.epsilon_t_prime_A = 1.0
        self.epsilon_t_prime_B = 1.0
        self.epsilon_t_prime_C = 1.0
        self.epsilon_t_prime_D = 1.0
        self.fcc_prime_A = 1.0
        self.fcc_prime_B = 1.0
        self.fcc_prime_C = 1.0
        self.fcc_prime_D = 1.0
        self.layers = 1
        self.alignment_depth = 1.0
        self.label_40.setVisible(False)
        self.lineEdit_7.setVisible(False)
        self.label_41.setVisible(False)
        self.radioButton.toggled.connect(self.toggle_visibility_rectangular_section)
        self.radioButton_2.toggled.connect(self.toggle_visibility_circular_section)
        self.pushButton.clicked.connect(self.data_retrival)
        self.pushButton.clicked.connect(self.calculations)
        self.pushButton.clicked.connect(self.point_A_prime)
        self.pushButton.clicked.connect(self.point_B)
        self.pushButton.clicked.connect(self.point_C)
        self.pushButton.clicked.connect(self.point_D)
        self.pushButton.clicked.connect(self.table_population)
        self.pushButton.clicked.connect(self.table_2_population)
        self.pushButton.clicked.connect(self.update_plot)
        self.x_input.editingFinished.connect(self.add_point)
        self.y_input.editingFinished.connect(self.add_point)
        self.pushButton_2.clicked.connect(self.generate_report)  # Connect button to function
        self.pushButton_3.clicked.connect(self.draw_cross_section)  # Add a new button to draw the cross-section
        self.canvas.mpl_connect('motion_notify_event', self.show_cursor_coordinates)

    def show_message(self):
        msg = QMessageBox(self)
        msg.setIcon(QMessageBox.Icon.Information)
        msg.setText('Name: Yousef Ahmed Mohamed Sedik\n\nID: 2001190\n\nFRP App and Report Assignment.')
        msg.setWindowTitle('About the developer.')
        msg.exec()

    def setup_canvas(self):
        layout = QVBoxLayout(self.frame_17)  # Ensure frame_17 is the correct widget
        layout.addWidget(self.canvas)  # Add the canvas to the layout
        self.create_canvas()  # Ensure this sets up the canvas

    def toggle_visibility_rectangular_section(self):
        is_checked = self.radioButton.isChecked()
        self.label.setVisible(is_checked)
        self.lineEdit.setVisible(is_checked)
        self.label_2.setVisible(is_checked)
        self.lineEdit_2.setVisible(is_checked)
        self.label_23.setVisible(is_checked)
        self.label_23.setVisible(is_checked)

    def toggle_visibility_circular_section(self):
        is_checked = self.radioButton_2.isChecked()
        self.label_40.setVisible(is_checked)
        self.lineEdit_7.setVisible(is_checked)
        self.label_41.setVisible(is_checked)

    def data_retrival(self):
        try:
            self.h = float(self.lineEdit.text()) if self.lineEdit.text() else 0
            self.b = float(self.lineEdit_2.text()) if self.lineEdit_2.text() else 0
            self.diameter = float(self.lineEdit_7.text()) if self.lineEdit_7.text() else 0
            self.r_c = float(self.lineEdit_3.text()) if self.lineEdit_3.text() else 0
            self.number_of_main_rebar = self.spinBox_4.value()
            self.dia_of_main_rebar = int(self.lineEdit_4.text()) if self.lineEdit_4.text() else 0
            self.dia_of_sec_rebar = int(self.lineEdit_5.text()) if self.lineEdit_5.text() else 0
            self.dia_of_side_rebar = int(self.lineEdit_6.text()) if self.lineEdit_6.text() else 0
            self.number_of_sec_rebar = self.spinBox_2.value()
            self.number_of_side_rebar = self.spinBox_5.value()
            self.fc_prime = float(self.lineEdit_8.text()) if self.lineEdit_8.text() else 0
            self.f_y = float(self.lineEdit_9.text()) if self.lineEdit_9.text() else 0
            self.E_s = float(self.lineEdit_10.text()) if self.lineEdit_10.text() else 200000.0
            self.f_fu_manu = float(self.lineEdit_11.text()) if self.lineEdit_11.text() else 0
            self.eps_fu_manu = float(self.lineEdit_12.text()) if self.lineEdit_12.text() else 0
            self.Ef = float(self.lineEdit_13.text()) if self.lineEdit_13.text() else 0
            self.t_f = float(self.lineEdit_14.text()) if self.lineEdit_14.text() else 0
            self.n_plies = self.spinBox.value()
            self.psi_f = float(self.lineEdit_19.text()) if self.lineEdit_19.text() else 0
            self.phi = float(self.lineEdit_21.text()) if self.lineEdit_21.text() else 0
            self.k_e = float(self.lineEdit_24.text()) if self.lineEdit_24.text() else 1
            self.cover = float(self.lineEdit_15.text()) if self.lineEdit_15.text() else 0
            exposure = self.comboBox.currentIndex()
            frp_mat = self.comboBox_2.currentIndex()
            CE_values = {(0, 0): 0.95, (0, 1): 0.85, (0, 2): 0.75, (1, 0): 0.85, (1, 1): 0.75, (1, 2): 0.65, (2, 0): 0.85, (2, 1): 0.7, (2, 2): 0.5}
            self.CE = CE_values.get((exposure, frp_mat))
            self.lineEdit_23.setPlaceholderText(str(self.CE))
        except Exception as e:
            print(f'Error, {e}')

    def calculations(self):
        try:
            self.A_g = self.h * self.b
            self.area_of_main_rebar = round(self.number_of_main_rebar * self.dia_of_main_rebar ** 2 / 4 * math.pi, 2)
            self.area_of_sec_rebar = round(self.number_of_sec_rebar * self.dia_of_sec_rebar ** 2 / 4 * math.pi, 2)
            self.area_of_side_rebar = round(self.number_of_side_rebar * self.dia_of_side_rebar ** 2 / 4 * math.pi, 2)
            self.A_s = self.area_of_main_rebar + self.area_of_sec_rebar + self.area_of_side_rebar
            self.f_fu = self.CE * self.f_fu_manu
            self.epsilon_fu = self.CE * self.eps_fu_manu
            self.epsilon_fe_POINT_A_prime = self.k_e * self.epsilon_fu
            self.epsilon_fe_POINT_B = min(0.004, self.k_e * self.epsilon_fu)
            self.epsilon_fe_POINT_C = min(0.004, self.k_e * self.epsilon_fu)
            self.epsilon_fe_POINT_D = min(0.004, self.k_e * self.epsilon_fu)
            self.diagonal = math.sqrt(self.h ** 2 + self.b ** 2)
            self.f_l_POINT_A_prime = 2 * self.n_plies * self.t_f * self.Ef * self.epsilon_fe_POINT_A_prime * self.psi_f / self.diagonal
            self.f_l_POINT_B = 2 * self.n_plies * self.t_f * self.Ef * self.epsilon_fe_POINT_B * self.psi_f / self.diagonal
            self.f_l_POINT_C = 2 * self.n_plies * self.t_f * self.Ef * self.epsilon_fe_POINT_C * self.psi_f / self.diagonal
            self.f_l_POINT_D = 2 * self.n_plies * self.t_f * self.Ef * self.epsilon_fe_POINT_D * self.psi_f / self.diagonal
            self.A_c = self.A_g - self.A_s
            self.rho_g = self.A_s / self.A_g
            self.A_e = self.A_g - ((self.h - 2 * self.r_c) ** 2 + (self.b - 2 * self.r_c) ** 2) / 3 - self.rho_g * self.A_g
            self.k_a = self.A_e / self.A_c * (self.b / self.h) ** 2
            self.k_b = self.A_e / self.A_c * (self.h / self.b) ** 0.5
            self.fcc_prime_A = self.fc_prime + 3.3 * self.psi_f * self.f_l_POINT_A_prime
            self.fcc_prime_B = self.fc_prime + 3.3 * self.psi_f * self.f_l_POINT_B
            self.fcc_prime_C = self.fc_prime + 3.3 * self.psi_f * self.f_l_POINT_C
            self.fcc_prime_D = self.fc_prime + 3.3 * self.psi_f * self.f_l_POINT_D
            self.epsilon_c_prime = 0.002
            self.epsilon_ccu_POINT_B = min(self.epsilon_c_prime * (1.5 + 12 * self.k_b * (self.f_l_POINT_B / self.fc_prime) * (self.epsilon_fe_POINT_B / self.epsilon_c_prime) ** 0.45), 0.01)
            self.epsilon_ccu_POINT_C = min(self.epsilon_c_prime * (1.5 + 12 * self.k_b * (self.f_l_POINT_C / self.fc_prime) * (self.epsilon_fe_POINT_C / self.epsilon_c_prime) ** 0.45), 0.01)
            self.epsilon_ccu_POINT_D = min(self.epsilon_c_prime * (1.5 + 12 * self.k_b * (self.f_l_POINT_D / self.fc_prime) * (self.epsilon_fe_POINT_D / self.epsilon_c_prime) ** 0.45), 0.01)
            self.d = self.h - self.cover
            self.c_point_C = self.d
            self.epsilon_f_y = self.f_y / self.E_s
            self.E_c = 4700 * np.sqrt(self.fc_prime)
            self.E_2_A = (self.fcc_prime_A - self.fc_prime) / self.epsilon_ccu_POINT_D
            self.E_2_B = (self.fcc_prime_B - self.fc_prime) / self.epsilon_ccu_POINT_B
            self.E_2_C = (self.fcc_prime_C - self.fc_prime) / self.epsilon_ccu_POINT_C
            self.E_2_D = (self.fcc_prime_D - self.fc_prime) / self.epsilon_ccu_POINT_D
            self.epsilon_t_prime_A = 2 * self.fc_prime / (self.E_c - self.E_2_A)
            self.epsilon_t_prime_B = 2 * self.fc_prime / (self.E_c - self.E_2_B)
            self.epsilon_t_prime_C = 2 * self.fc_prime / (self.E_c - self.E_2_C)
            self.epsilon_t_prime_D = 2 * self.fc_prime / (self.E_c - self.E_2_D)
            self.layers = int(2 + self.number_of_side_rebar / 2)
            self.bars_spacing = (self.h - 2 * self.cover) / self.layers
            self.alignment_depth = self.h - 2 * self.cover
            d_i = [self.alignment_depth / 2 - self.bars_spacing * layer for layer in range(self.layers + 1)]
        except ZeroDivisionError:
            print('Can\'t divide by zero')

    def point_A_prime(self):
        try:
            Pn_A_prime = point_A_prime_function(self.phi, self.fcc_prime_A, self.A_g, self.A_s, self.f_y)
            Mn_A_prime = 0
            return (Pn_A_prime, Mn_A_prime)
        except Exception:
            pass

    def point_B(self):
        try:
            c = self.d
            var_A = A_function(self.b, self.E_c, self.E_2_B, self.fc_prime, self.epsilon_ccu_POINT_B, c)
            var_B = B_function(self.b, self.E_c, self.E_2_B, self.epsilon_ccu_POINT_B, c)
            var_C = C_function(self.b, self.fc_prime)
            var_D = D_function(self.b, c, self.fc_prime, self.E_2_B, self.epsilon_ccu_POINT_B)
            var_E = E_function(self.b, self.E_c, self.E_2_B, self.fc_prime, self.epsilon_ccu_POINT_B, c)
            var_F = F_function(self.b, self.h, c, self.E_c, self.E_2_B, self.fc_prime, self.epsilon_ccu_POINT_B)
            var_G = G_function(self.b, self.h, c, self.E_c, self.E_2_B, self.fc_prime, self.epsilon_ccu_POINT_B)
            var_H = H_function(self.b, self.h, c, self.fc_prime)
            var_I = I_function(self.b, self.h, c, self.E_2_B, self.fc_prime, self.epsilon_ccu_POINT_B)
            yt = c * self.epsilon_t_prime_B / self.epsilon_ccu_POINT_B
            d_si = [(self.d - self.cover) / 2 - self.bars_spacing * layer for layer in range(self.layers + 1)]
            f_si = [steel_strain(self.d, self.cover, self.bars_spacing, layer, self.epsilon_ccu_POINT_B, self.E_s, self.f_y) * self.E_s for layer in range(self.layers + 1)]
            A_s_top = [self.dia_of_sec_rebar ** 2 / 4 * 3.14 * self.number_of_sec_rebar]
            A_s_side = [self.dia_of_side_rebar ** 2 / 4 * 3.14 * 2 for i in range(1, self.layers)]
            A_s_bottom = [self.dia_of_main_rebar ** 2 / 4 * 3.14 * self.number_of_main_rebar]
            A_si = A_s_top + A_s_side + A_s_bottom
            side_bars_forces = sum((A_si[i] * f_si[i] * d_si[i] for i in range(self.layers)))
            A_s = sum(A_si)
            Pn_point_B = Pn_function(var_A, var_B, var_C, var_D, A_s, self.f_y, yt)
            Mn_point_B = Mn_function(var_E, var_F, var_G, var_H, var_I, side_bars_forces, yt)
            return (Pn_point_B, Mn_point_B)
        except Exception:
            pass

    def point_C(self):
        try:
            c = self.d * (self.epsilon_ccu_POINT_C / (self.epsilon_f_y + self.epsilon_ccu_POINT_C))
            var_A = A_function(self.b, self.E_c, self.E_2_C, self.fc_prime, self.epsilon_ccu_POINT_C, c)
            var_B = B_function(self.b, self.E_c, self.E_2_C, self.epsilon_ccu_POINT_C, c)
            var_C = C_function(self.b, self.fc_prime)
            var_D = D_function(self.b, c, self.fc_prime, self.E_2_C, self.epsilon_ccu_POINT_C)
            var_E = E_function(self.b, self.E_c, self.E_2_C, self.fc_prime, self.epsilon_ccu_POINT_C, c)
            var_F = F_function(self.b, self.h, c, self.E_c, self.E_2_C, self.fc_prime, self.epsilon_ccu_POINT_C)
            var_G = G_function(self.b, self.h, c, self.E_c, self.E_2_C, self.fc_prime, self.epsilon_ccu_POINT_C)
            var_H = H_function(self.b, self.h, c, self.fc_prime)
            var_I = I_function(self.b, self.h, c, self.E_2_C, self.fc_prime, self.epsilon_ccu_POINT_C)
            yt = c * self.epsilon_t_prime_C / self.epsilon_ccu_POINT_C
            d_si = [(self.d - self.cover) / 2 - self.bars_spacing * layer for layer in range(self.layers + 1)]
            f_si = [steel_strain(self.d, self.cover, self.bars_spacing, layer, self.epsilon_ccu_POINT_C, self.E_s, self.f_y) * self.E_s for layer in range(self.layers + 1)]
            A_s_top = [self.dia_of_sec_rebar ** 2 / 4 * 3.14 * self.number_of_sec_rebar]
            A_s_side = [self.dia_of_side_rebar ** 2 / 4 * 3.14 * 2 for i in range(1, self.layers)]
            A_s_bottom = [self.dia_of_main_rebar ** 2 / 4 * 3.14 * self.number_of_main_rebar]
            A_si = A_s_top + A_s_side + A_s_bottom
            side_bars_forces = sum((A_si[i] * f_si[i] * d_si[i] for i in range(self.layers)))
            A_s = sum(A_si)
            Pn_point_C = Pn_function(var_A, var_B, var_C, var_D, A_s, self.f_y, yt)
            Mn_point_C = Mn_function(var_E, var_F, var_G, var_H, var_I, side_bars_forces, yt)
            return (Pn_point_C, Mn_point_C)
        except Exception:
            pass

    def point_D(self):
        try:
            c = 0.1 * self.d / 0.8
            var_E = E_function(self.b, self.E_c, self.E_2_D, self.fc_prime, self.epsilon_ccu_POINT_D, c)
            var_F = F_function(self.b, self.h, c, self.E_c, self.E_2_D, self.fc_prime, self.epsilon_ccu_POINT_D)
            var_G = G_function(self.b, self.h, c, self.E_c, self.E_2_D, self.fc_prime, self.epsilon_ccu_POINT_D)
            var_H = H_function(self.b, self.h, c, self.fc_prime)
            var_I = I_function(self.b, self.h, c, self.E_2_D, self.fc_prime, self.epsilon_ccu_POINT_D)
            yt = c * self.epsilon_t_prime_D / self.epsilon_ccu_POINT_D
            d_si = [(self.d - self.cover) / 2 - self.bars_spacing * layer for layer in range(self.layers + 1)]
            f_si = [steel_strain(self.d, self.cover, self.bars_spacing, layer, self.epsilon_ccu_POINT_C, self.E_s, self.f_y) * self.E_s for layer in range(self.layers + 1)]
            A_s_top = [self.dia_of_sec_rebar ** 2 / 4 * 3.14 * self.number_of_sec_rebar]
            A_s_side = [self.dia_of_side_rebar ** 2 / 4 * 3.14 * 2 for i in range(1, self.layers)]
            A_s_bottom = [self.dia_of_main_rebar ** 2 / 4 * 3.14 * self.number_of_main_rebar]
            A_si = A_s_top + A_s_side + A_s_bottom
            side_bars_forces = sum((A_si[i] * f_si[i] * d_si[i] for i in range(self.layers)))
            Pn_point_C = 0
            Mn_point_C = Mn_function(var_E, var_F, var_G, var_H, var_I, side_bars_forces, yt)
            return (Pn_point_C, Mn_point_C)
        except Exception:
            pass

    def table_population(self):
        try:
            self.tableWidget.setItem(0, 0, QTableWidgetItem(f'{(self.h if self.radioButton.isChecked() else 0)} mm'))
            self.tableWidget.setItem(0, 1, QTableWidgetItem('Height of the cross-section.'))
            self.tableWidget.setItem(1, 0, QTableWidgetItem(f'{(self.b if self.radioButton.isChecked() else 0)} mm'))
            self.tableWidget.setItem(1, 1, QTableWidgetItem('Width of the cross-section.'))
            self.tableWidget.setItem(2, 0, QTableWidgetItem(f'{(self.diameter if self.radioButton_2.isChecked() else 0)} mm'))
            self.tableWidget.setItem(2, 1, QTableWidgetItem('Diameter of the cross-section.'))
            self.tableWidget.setItem(3, 0, QTableWidgetItem(f'{self.r_c} mm'))
            self.tableWidget.setItem(3, 1, QTableWidgetItem('Radius of the corners.'))
            self.tableWidget.setItem(4, 0, QTableWidgetItem(f'{self.A_g} mm²'))
            self.tableWidget.setItem(4, 1, QTableWidgetItem(f'{self.b} x {self.h}.'))
            self.tableWidget.setItem(5, 0, QTableWidgetItem(f'{self.area_of_main_rebar} mm²'))
            self.tableWidget.setItem(5, 1, QTableWidgetItem(f'{self.number_of_main_rebar} x ({self.dia_of_main_rebar}²π/4)'))
            self.tableWidget.setItem(6, 0, QTableWidgetItem(f'{self.area_of_sec_rebar} mm²'))
            self.tableWidget.setItem(6, 1, QTableWidgetItem(f'{self.number_of_sec_rebar} x ({self.dia_of_sec_rebar}²π/4)'))
            self.tableWidget.setItem(7, 0, QTableWidgetItem(f'{self.area_of_side_rebar} mm²'))
            self.tableWidget.setItem(7, 1, QTableWidgetItem(f'{self.number_of_side_rebar} x ({self.dia_of_side_rebar}²π/4)'))
            self.tableWidget.setItem(8, 0, QTableWidgetItem(f'{self.fc_prime} MPa'))
            self.tableWidget.setItem(8, 1, QTableWidgetItem('User defined.'))
            self.tableWidget.setItem(9, 0, QTableWidgetItem(f'{self.f_fu_manu} MPa'))
            self.tableWidget.setItem(9, 1, QTableWidgetItem('User defined.'))
            self.tableWidget.setItem(10, 0, QTableWidgetItem(f'{self.eps_fu_manu} mm/mm'))
            self.tableWidget.setItem(10, 1, QTableWidgetItem('User defined.'))
            self.tableWidget.setItem(11, 0, QTableWidgetItem(f'{self.Ef} MPa'))
            self.tableWidget.setItem(11, 1, QTableWidgetItem('User defined.'))
            self.tableWidget.setItem(12, 0, QTableWidgetItem(f'{self.t_f} mm'))
            self.tableWidget.setItem(12, 1, QTableWidgetItem('User defined.'))
            self.tableWidget.setItem(13, 0, QTableWidgetItem(f'{self.n_plies}'))
            self.tableWidget.setItem(13, 1, QTableWidgetItem('User defined.'))
            self.tableWidget.setItem(14, 0, QTableWidgetItem(f'{self.CE}'))
            self.tableWidget.setItem(14, 1, QTableWidgetItem(f'Environment reduction factor according to ACI440.2R-2017\nExposure type: {self.comboBox.currentText()} and FRP material: {self.comboBox_2.currentText()}0'))
            self.tableWidget.setItem(15, 0, QTableWidgetItem(f'{self.psi_f}'))
            self.tableWidget.setItem(15, 1, QTableWidgetItem('Further reduction factor according to ACI440.2R-2017'))
            self.tableWidget.setItem(16, 0, QTableWidgetItem(f'{self.phi}'))
            self.tableWidget.setItem(16, 1, QTableWidgetItem('Strength reduction factor according to ACI440.2R-2017'))
            self.tableWidget.setItem(17, 0, QTableWidgetItem(f'{self.f_fu} MPa'))
            self.tableWidget.setItem(17, 1, QTableWidgetItem(f'{self.CE} x {self.f_fu_manu}'))
            self.tableWidget.setItem(18, 0, QTableWidgetItem(f'{self.epsilon_fu} mm/mm'))
            self.tableWidget.setItem(18, 1, QTableWidgetItem(f'{self.CE} x {self.eps_fu_manu}'))
            self.tableWidget.setItem(19, 0, QTableWidgetItem(f'{self.k_e}'))
            self.tableWidget.setItem(19, 1, QTableWidgetItem('FRP strain efficiency factor. Average 0.586'))
            self.tableWidget.setItem(20, 0, QTableWidgetItem(f'{self.point_A_prime()[0] / 1000.0}0 kN'))
            self.tableWidget.setItem(21, 0, QTableWidgetItem(f'{self.point_A_prime()[1] / 1000000.0} kN-m'))
            self.tableWidget.setItem(22, 0, QTableWidgetItem(f'{self.point_B()[0] / 1000.0}0 kN'))
            self.tableWidget.setItem(23, 0, QTableWidgetItem(f'{self.point_B()[1] / 1000000.0} kN-m'))
            self.tableWidget.setItem(24, 0, QTableWidgetItem(f'{self.point_C()[0] / 1000.0}0 kN'))
            self.tableWidget.setItem(25, 0, QTableWidgetItem(f'{self.point_C()[1] / 1000000.0} kN-m'))
            self.tableWidget.setItem(26, 0, QTableWidgetItem(f'{self.point_D()[0] / 1000.0}0 kN'))
            self.tableWidget.setItem(27, 0, QTableWidgetItem(f'{self.point_D()[1] / 1000000.0} kN-m'))
            self.tableWidget.resizeColumnsToContents()
            self.tableWidget.resizeRowsToContents()
            self.tableWidget.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
            self.tableWidget.setAlternatingRowColors(True)
        except Exception:
            pass

    def table_2_population(self):
        try:
            self.tableWidget_2.setItem(0, 0, QTableWidgetItem(str(round(self.point_A_prime()[1] / 1000000.0, 2))))
            self.tableWidget_2.setItem(0, 1, QTableWidgetItem(str(round(self.point_A_prime()[0] / 1000.0, 2))))
            self.tableWidget_2.setItem(1, 0, QTableWidgetItem(str(round(0.65 * self.point_B()[1] / 1000000.0, 2))))
            self.tableWidget_2.setItem(1, 1, QTableWidgetItem(str(round(0.65 * self.point_B()[0] / 1000.0, 2))))
            self.tableWidget_2.setItem(2, 0, QTableWidgetItem(str(round(self.point_C()[1] / 1000000.0, 2))))
            self.tableWidget_2.setItem(2, 1, QTableWidgetItem(str(round(self.point_C()[0] / 1000.0, 2))))
            self.tableWidget_2.setItem(3, 0, QTableWidgetItem(str(round(self.point_D()[1] / 1000000.0, 2))))
            self.tableWidget_2.setItem(3, 1, QTableWidgetItem(str(round(self.point_D()[0] / 1000.0, 2))))
            self.tableWidget_2.setItem(4, 0, QTableWidgetItem(str(round(float(self.x_input.text()), 2))))
            self.tableWidget_2.setItem(4, 1, QTableWidgetItem(str(round(float(self.y_input.text()), 2))))
            self.tableWidget_2.resizeColumnsToContents()
            self.tableWidget_2.resizeRowsToContents()
            self.tableWidget_2.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
            self.tableWidget_2.setAlternatingRowColors(True)
        except Exception:
            pass

    def create_canvas(self):
        try:
            self.ax.clear()  # Clear previous plots
            self.ax.set_title('FRP ID')
            self.ax.set_xlabel('Moment in kN-m')
            self.ax.set_ylabel('Normal in kN')
            self.ax.axhline(color='black', linestyle='-', linewidth=1)
            self.ax.axvline(color='black', linestyle='-', linewidth=1)
            self.canvas.draw()  # Redraw the canvas
        except Exception as e:
            print(f"Error in create_canvas: {e}")

    def add_point(self):
        try:
            x = float(self.x_input.text())
            y = float(self.y_input.text())

            # Append the new point to the list (do not clear)
            self.points.clear()
            self.points.append((x, y))

            # Plot the new point
            self.ax.clear()
            # self.create_canvas()
            self.ax.plot(x, y, 'ro', label='User-defined.' if len(self.points) == 1 else "")
            self.ax.legend()

            # Redraw the canvas
            self.canvas.draw()

        except ValueError:
            print("Invalid input")  # Debugging line

    def update_plot(self):
        try:
            self.ax.clear()  # Clear the previous plot
            self.create_canvas()  # Recreate canvas

            # Example points - you need to ensure these methods return valid data
            x = [self.point_A_prime()[1] / 1000000.0, 0.65 * self.point_B()[1] / 1000000.0, 0.65 * self.point_C()[1] / 1000000.0, 0.65 * self.point_D()[1] / 1000000.0]
            y = [self.point_A_prime()[0] / 1000.0, 0.65 * self.point_B()[0] / 1000.0, 0.65 * self.point_C()[0] / 1000.0, 0.65 * self.point_D()[0] / 1000.0]

            # Plot the main line
            self.ax.plot(x, y, label='FRP ID')

            # Plot user-defined points
            if self.points:
                user_x, user_y = zip(*self.points)
                self.ax.plot(user_x, user_y, 'ro', label='User-defined.')
                for (ux, uy) in zip(user_x, user_y):
                    self.ax.text(ux, uy, f'({ux:.1f}, {uy:.1f})', fontsize=11, color='red', ha='center', va='bottom')

            self.ax.legend()  # Add legend
            self.canvas.draw()  # Redraw the canvas

        except Exception as e:
            print(f"Error in update_plot: {e}")

    def show_cursor_coordinates(self, event):
        try:
            if event.inaxes == self.ax:
                self.cursor_label.setText(f"Cursor: ({event.xdata:.2f}, {event.ydata:.2f})")
        except Exception:
            pass

    def generate_report(self):
        try:
            # Set up the PDF document
            report_filename = "FRP_Report.pdf"
            c = canvas.Canvas(report_filename, pagesize=letter)
            width, height = letter  # Default page size

            # Title of the report
            c.setFont("Helvetica-Bold", 16)
            c.drawString(200, height - 40, "FRP-confined RC Column Report")

            # Add the FRP ID plot to the PDF (assuming you've saved the plot as an image)
            plot_image_path = "plot_image.png"
            self.canvas.print_figure(plot_image_path, bbox_inches='tight')

            # Add plot image to the PDF
            c.drawImage(plot_image_path, 50, height - 300, width=500, height=200)

            # Add data to the PDF
            c.setFont("Helvetica", 12)
            c.drawString(50, height - 320, "Column Details:")

            # Example of adding dynamic data (use your actual variables)
            c.drawString(50, height - 340, f"Column Height (h): {self.h} m")
            c.drawString(50, height - 360, f"Column Width (b): {self.b} m")
            c.drawString(50, height - 380, f"Main Rebar Diameter: {self.dia_of_main_rebar} mm")
            c.drawString(50, height - 400, f"Concrete Strength (f'c): {self.fc_prime} MPa")

            # User-defined points data (if available)
            if self.points:
                c.drawString(50, height - 420, "User-Defined Straining Actions:")
                y_offset = height - 440
                for i, (ux, uy) in enumerate(self.points):
                    c.drawString(50, y_offset - (i * 20), f"Point {i+1}: ({ux:.1f}, {uy:.1f}) kN-m, kN")
                    if y_offset - (i * 20) < 50:  # Check if we're getting close to the bottom of the page
                        c.showPage()  # Start a new page
                        y_offset = height - 40

            # Save the PDF
            c.save()
            print(f"Report saved as {report_filename}")
        except Exception as e:
            print(f"Error generating report: {e}")

    def draw_cross_section(self):
        # Column dimensions
        column_width = self.b  # Width of the column
        column_height = self.h  # Height of the column

        # Rebar details (you can change this to take values from user input)
        main_rebar_diameter = self.dia_of_main_rebar  # Diameter of the main rebar
        num_main_rebars = self.number_of_main_rebar  # Number of main rebars

        # Create column cross-section (rectangle)
        column = Polygon([(0, 0), (column_width, 0), (column_width, column_height), (0, column_height)])

        # Create main rebars (circles)
        rebar_circles = []
        spacing = column_width / (num_main_rebars + 1)  # Simple uniform spacing for rebars
        for i in range(1, num_main_rebars + 1):
            rebar_center = (spacing * i, column_height / 2)  # Placing rebars at the center of the column
            rebar = Point(rebar_center).buffer(main_rebar_diameter / 2)  # Create circle for rebar
            rebar_circles.append(rebar)

        # Create plot
        fig, ax = plt.subplots()
        ax.set_aspect('equal', 'box')
        ax.set_title('Column Cross-Section with Rebar')
        ax.set_xlabel('Width (m)')
        ax.set_ylabel('Height (m)')

        # Plot the column
        column_patch = MPolygon(np.array(column.exterior.xy).T, closed=True, edgecolor='black', facecolor='gray', alpha=0.3)
        ax.add_patch(column_patch)

        # Plot the rebars
        for rebar in rebar_circles:
            # We manually pass the radius of the rebar here
            rebar_patch = MCircle(rebar.centroid, main_rebar_diameter / 2, edgecolor='red', facecolor='red', alpha=0.6)
            ax.add_patch(rebar_patch)

        # Show the plot
        plt.show()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainApp()
    window.show()
    window.show_message()
    app.exec()