import sys
import math
import random
# The function 'uic' allow us to load the GUI - Graphic User Interface
from PyQt5 import uic
from PyQt5.QtCore import QDir

from PyQt5.QtWidgets import (
    QMainWindow,  # QMainWindow is our current window
    QApplication, # QApplication is used to load the application
    QFileDialog,
    QMessageBox
)
# Local functions
from utils import *

# Load the file interface.ui
qtCreatorFile = "interface.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

n = 1


class Modelo():

    def __init__(self):
        self.a_sup = 0
        self.a_inf = 0
        self.b_sup = 0
        self.b_inf = 0
        self.c_sup = 0
        self.c_inf = 0

    def validateA(self, value):
        return  True if self.a_inf <= value <= self.a_sup else False

    def validateB(self, value):
        return  True if self.b_inf <= value <= self.b_sup else False
    
    def validateC(self, value):
        return  True if self.c_inf <= value <= self.c_sup else False
        

class Application(QMainWindow, Ui_MainWindow):

    def __init__(self):
         # Initialize all the constructors
        QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        # Set a fixed size for the window using: setFixedSize(width, height) 
        QMainWindow.setFixedSize(self, 1109, 828)
        # Class variables
        self.x = []
        self.y = []

        """ Buttons's Events """
        self.btn_loadfile.clicked.connect(self.loadFiles)
        self.btn_start.clicked.connect(self.start)

    
    def loadFiles(self):
        file, _ = QFileDialog.getOpenFileName(self, 'Select a File', QDir.currentPath(), "Text Files (*.txt)")
        # If the file exists
        if file:
            self.fname, self.ext = getFileExtension(file)
            self.txt_file.setText(f'{self.fname}.{self.ext}')
            # Reading the content of the chosen file line by line
            lines = open(file, 'r', encoding='utf-8').readlines()
            # Strip the newline ‘\n’ character using strip() function 
            for line in lines:
                x = int(line.strip().split(',')[0])
                y = int(line.strip().split(',')[1])
                # Add each coordinate to their respective array
                self.x.append(x)
                self.y.append(y)

            # Maximum values of x and y 
            self.maximum_x = max(list(map( lambda number: abs(number), self.x )))
            self.maximum_y = max(list(map( lambda number: abs(number), self.y )))
            self.maximum = max(self.maximum_x, self.maximum_y)

            # Compute constraints for lineal equation
            self.l = Modelo()

            self.l.a_sup = 20
            self.l.a_inf = -20
            self.txt_las.setText(str(self.l.a_sup))
            self.txt_lai.setText(str(self.l.a_inf))

            self.l.b_sup = 2 * self.maximum_y
            self.l.b_inf = -2 * self.maximum_y
            self.txt_lbs.setText(str(self.l.b_sup)) 
            self.txt_lbi.setText(str(self.l.b_inf))

            # Compute constraints for cuadratic equation
            self.c = Modelo()

            self.c.a_sup = self.maximum 
            self.c.a_inf = -self.maximum 
            self.txt_cas.setText(str(self.c.a_sup))
            self.txt_cai.setText(str(self.c.a_inf))

            self.c.b_sup = 2*self.maximum 
            self.c.b_inf = -2*self.maximum 
            self.txt_cbs.setText(str(self.c.b_sup)) 
            self.txt_cbi.setText(str(self.c.b_inf))

            self.c.c_sup = 2 * self.maximum_y
            self.c.c_inf = -2 * self.maximum_y
            self.txt_ccs.setText(str(self.c.c_sup)) 
            self.txt_cci.setText(str(self.c.c_inf)) 

            # Compute constraints for gaussian equation
            self.g = Modelo()

            self.g.a_sup = 2 * self.maximum_y
            self.g.a_inf = -2 * self.maximum_y
            self.txt_gas.setText(str(self.g.a_sup)) 
            self.txt_gai.setText(str(self.g.a_inf))
             
            self.g.b_sup = self.maximum_x
            self.txt_gbs.setText(str(self.g.b_sup)) 
            self.txt_gbi.setText(str(self.g.b_inf)) 

            self.g.c_sup = self.maximum_x
            self.g.c_inf = -self.maximum_x
            self.txt_gcs.setText(str(self.g.c_sup)) 
            self.txt_gci.setText(str(self.g.c_inf))  

        else: 
            self.createMessageBox('Error!', 'The selected file cannot be opened')

    def start(self):

        # Get the number of iterations
        iterations  = int(self.txt_iteration.toPlainText())
        populations = int(self.txt_n.toPlainText())
        individuals = int(self.txt_m.toPlainText())

        mjla = math.ceil(math.log(2*self.l.a_sup*(10**n))/math.log(2))
        mjc  = math.ceil(math.log(2*self.l.b_sup*(10**n))/math.log(2))
        mjca = math.ceil(math.log(2*self.c.a_sup*(10**n))/math.log(2))
        mjcb = math.ceil(math.log(2*self.c.b_sup*(10**n))/math.log(2))
        mjgb = math.ceil(math.log(self.g.b_sup*(10**n))/math.log(2))
        mjgc = math.ceil(math.log(2*self.g.c_sup*(10**n))/math.log(2))

        cromosomas_lineal     = mjla + mjc        
        cromosomas_cuadratica = mjca + mjcb + mjc 
        cromosomas_gaussiana  = mjc + mjgb + mjgc 
        
        index_1 = 1
        best_vector_1 = computeLinear( populations, individuals, self.l, cromosomas_lineal, mjla, mjc, self.x, self.y )

        print('\t\t FUNCIÓN LINEAL \t\t\n')
        print('\nITERACION \t\t A \t\t\t B \t\t\t Z\n')
        print(f'  [{index_1}]  \t {best_vector_1[0]} \t {best_vector_1[1]} \t {best_vector_1[2]}\n')
        
        for i in range(iterations-1):
            current_vector = computeLinear( populations, individuals, self.l, cromosomas_lineal, mjla, mjc, self.x, self.y )
            print('\nITERACION \t\t A \t\t\t B \t\t\t Z')
            print(f'  [{i+2}]  \t {current_vector[0]} \t {current_vector[1]} \t {current_vector[2]}')
            if current_vector[-1] > best_vector_1[-1]:
                best_vector_1 = current_vector
                index_1 = i + 2

        index_2 = 1
        best_vector_2 = computeCG( populations, individuals, self.c, cromosomas_cuadratica, mjca, mjcb, mjc, self.x, self.y, True )
        print('\n\t\t FUNCIÓN CUADRÁTICA \t\t\n')
        print('\nITERACION \t\t A \t\t\t B \t\t\t C \t\t\t Z')
        print(f'  [{index_2}]  \t {best_vector_2[0]} \t {best_vector_2[1]} \t {best_vector_2[2]} \t {best_vector_2[3]}\n')
        
        for i in range(iterations-1):
            current_vector = computeCG( populations, individuals, self.c, cromosomas_cuadratica, mjca, mjcb, mjc, self.x, self.y, True )
            print('\nITERACION \t\t A \t\t\t B \t\t\t C \t\t\t Z')
            print(f'  [{i+2}]  \t {current_vector[0]} \t {current_vector[1]} \t {current_vector[2]} \t {current_vector[3]}')
            if current_vector[-1] > best_vector_2[-1]:
                best_vector_2 = current_vector
                index_2 = i + 2

        index_3 = 1
        best_vector_3 = computeCG( populations, individuals, self.g, cromosomas_gaussiana, mjc, mjgb, mjgc, self.x, self.y, False )
        print('\n\t\t FUNCIÓN GAUSSIANA \t\t\n')
        print('\nITERACION \t\t A \t\t\t B \t\t\t C \t\t\t Z')
        print(f'  [{index_3}]  \t {best_vector_3[0]} \t {best_vector_3[1]} \t {best_vector_3[2]} \t {best_vector_3[3]}\n')
        
        for i in range(iterations-1):
            current_vector = computeCG( populations, individuals, self.g, cromosomas_gaussiana, mjc, mjgb, mjgc, self.x, self.y, False )
            print('\nITERACION \t\t A \t\t\t B \t\t\t C \t\t\t Z')
            print(f'  [{i+2}]  \t {current_vector[0]} \t {current_vector[1]} \t {current_vector[2]} \t {current_vector[3]}')
            if current_vector[-1] > best_vector_3[-1]:
                best_vector_3 = current_vector
                index_3 = i + 2
        
        print(f'\nMejor Vector Función Lineal: {best_vector_1} en la iteración: {index_1}')
        print(f'\nMejor Vector Función Cuadrática: {best_vector_2} en la iteración: {index_2}')
        print(f'\nMejor Vector Función Gaussiana: {best_vector_3} en la iteración: {index_3}')



    def createMessageBox(self, title, maintext):
        box = QMessageBox()
        box.setIcon(QMessageBox.Information)
        box.setWindowTitle(title)
        box.setText(maintext)
        # Mostramos la caja de texto
        box.exec_()

# Cuando ejecutemos la app todo lo que esté aquí definido se ejecutará
if __name__ == '__main__':
    # Iniciamos la aplicación
    app = QApplication(sys.argv)
    # Inicializamos nuestra clase
    GUI = Application()
    # Mostramos nuestra aplicación
    GUI.show()
    # Cerramos nuestra aplicación
    sys.exit(app.exec_())
