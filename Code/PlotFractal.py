#!/usr/bin/env sage -python

######################################################################
# Autor: Andres Herrera Poyatos
# Universidad de Granada, Marzo, 2015
# Fractales de Newton
#######################################################################

# Ejecucion:
#  sage -python PlotFractal.py <polinomio> <densidad> <iteraciones Newton-Raphson> <tolerancia Newton-Raphson> <nombre de la imagen de salida>
# Densidad indica la dimensión de la matriz de puntos sobre los que se aplica Newton-Raphson.
# Ejemplo:
#  sage -python PlotFractal.py x**3-1 200 40 0.001 prueba.png

from sage.all import * # Se importa toda la funcionalidad de sage
import numpy as np     # Se importa la libreria numpy
import sys             # sys.argv
import time            # time para medir el tiempo

#-----------------------------------------------------------------------------------------------#

# Metodo de Newton-Raphson.
# Realiza n iteraciones del metodo y devuelve
# el resultado obtenido que se espera que sea una raiz de f.
# Parametros:
#  - f : Funcion sobre la que se aplica el metodo.
#  - f_derivada : Derivada de la funcion f.
#  - x_0 : Numero complejo tomado como aproximacion inicial.
#          Segun este numero se convergera a una raiz u otras.
#  - n : Numero de iteraciones a realizar.
#  - tolerancia : Si |x_n+1 - x_n| < tolerancia se para la ejecucion.
#      Se presupone que ya ha convergido.
#  Return:
#    Raiz obtenida (es una aproximacion).
def metodoNewtonRaphson(f , f_derivada, x_0, n, tolerancia):
    sol = x_0
    for i in xrange(0 , n ):
        derivada = N(f_derivada(sol))
        if derivada == 0:
            return sol, i
        prev = sol
        sol = sol - N(f(sol)) / derivada
        if (prev - sol).norm() < tolerancia:
            return sol, i
    
    return sol, n

#-----------------------------------------------------------------------------------------------#

# PlotFractal.
# Funcion que dibuja el fractal resultado de estudiar las raices a
# las que convergen mediante el metodo de Newton-Raphson los elementos
# del conjunto de los numeros complejos.
# Parametros:
# - expresion : String con la funcion de la que se dibuja el fractal asociado.
# - densidad : Dimension del grid sobre el que se realizaran las medidas. A mayor
#              densidad la imagen final sera de mayor calidad.
# - max_iteraciones : Numero maximo de iteraciones que realiza el metodo de Newton-Raphson.
# - tolerancia :  Si |x_n+1 - x_n| < tolerancia se para la ejecucion del metodo de Newton-Raphson.
#      Se presupone que ya ha convergido.
# - imagen : Nombre de la imagen donde se guarda el fractal.
def PlotFractal(expresion, densidad, max_iteraciones, tolerancia, imagen):
    # Se define la variable y se obtiene la funcion f
    x = var('x')
    f = eval(expresion)

    # Se calcula la derivada de f y sus raices
    f_derivada = diff(f, x)
    raices = [N(raiz.rhs()) for raiz in solve(f, x)]
    

    # Se obtiene un grid del cuadrado [-1,1]x[-1,1]
    n = np.linspace(-1, 1, num=densidad)
    [X,Y] = np.meshgrid(n,n)

    # Lista de listas donde se almacenaran los valores del fractal
    fractal = [[-1 for j in xrange(0,densidad) ] for i in xrange(0,densidad)]
    
    # Se calculan las raices a las que convergen los numeros complejos X[i,j]+ Y[i,j]i
    for i in xrange(0, densidad):
        for j in xrange(0,densidad):
 
            # Se calcula la raiz a la que converge X[i,j]+ Y[i,j]i y el numero de iteraciones requerido
            [raiz, k] = metodoNewtonRaphson(f, f_derivada, CDF(X[i,j], Y[i,j]), max_iteraciones, tolerancia)

            # Se obtiene el indice de la raiz a la que ha convergido
            min_indice = 0; min_valor = (raiz-raices[0]).norm()
            for d in range(1,len(raices)):
                valor = (raiz-raices[d]).norm()
                if valor < min_valor:
                    min_valor = valor; min_indice = d

            # Se asigna el valor correspondiente al fractal
            fractal[i][j] = min_indice + k*0.025

    F = matrix(fractal)
    p = F.plot(cmap = 'brg_r')
    p.save(imagen, title="Fractal de Newton para " + expresion)

#-----------------------------------------------------------------------------------------------#

######################## MAIN ##########################

# Comprobacion de los argumentos
if len(sys.argv) != 6:
    print("Sintaxis: sage -python PlotFractal <polinomio> <densidad> <iteraciones Newton-Raphson> <tolerancia Newton-Raphson> <nombre de la imagen de salida>")
    sys.exit()

# Llamada a la funcion
start_time = time.time()
PlotFractal(expresion=sys.argv[1], densidad=int(sys.argv[2]), max_iteraciones=int(sys.argv[3]), tolerancia=float(sys.argv[4]), imagen=sys.argv[5])
print("--- %f segundos ---" % (time.time() - start_time) )
print("Se ha obtenido el fractal " + sys.argv[5])

