import math
import random
from typing import Text

""" This function returns back the extension and the file name of the given path """
def getFileExtension(file_path):
    snippets = file_path.split('/')[-1].split('.')
    fname = snippets[0]
    ext = snippets[1]
    return (fname, ext)

""" This function takes in an 'int' and returns back the number as binary representation """
def generateBits( bits ):
    number = random.getrandbits(bits)
    binary = bin(number).split('b')[1]
    list_binary = list(binary)

    while len(list_binary) < bits:
        list_binary.insert(0, '0')
    
    return "".join(list_binary)

def computePhenotype(lim_sup, lim_inf, substring, bits):
    # print(f'SUBSTRING: {substring}')
    number = int(substring, 2)
    return lim_inf + number * ( (lim_sup - lim_inf) / ((2**bits) -1) )

""" Esta función genera los fenotipos de la cadena de bist dada """
def generatePhenotypesLineal(bits, linear, mja, mjb):

    phenotype_a = 0.0
    phenotype_b = 0.0

    phenotype_a = computePhenotype(linear.a_sup, linear.a_inf, bits[:mja], mja)
    phenotype_b = computePhenotype(linear.b_sup, linear.b_inf, bits[mja:], mjb)
    val = linear.validateA(phenotype_a) and linear.validateB(phenotype_b)

    return val, phenotype_a, phenotype_b

def generatePhenotypesCG(bits, gc, mja, mjb, mjc):

    phenotype_a = 0.0 
    phenotype_b = 0.0 
    phenotype_c = 0.0 

    phenotype_a = computePhenotype(gc.a_sup, gc.a_inf, bits[:mja], mja)
    phenotype_b = computePhenotype(gc.b_sup, gc.b_inf, bits[mja:mja+mjb], mjb)
    phenotype_c = computePhenotype(gc.c_sup, gc.c_inf, bits[mja+mjb:], mjc)
    val = gc.validateA(phenotype_a) and gc.validateB(phenotype_b) and gc.validateC(phenotype_c)

    return val, phenotype_a, phenotype_b, phenotype_c 

def computeLinear( populations, individuals, linear, num_cromosmas, mja, mjb, x, y ):

    linear_vectors = dict()
    best_vector = []
    
    for i in range(populations):
        competent_vectors = []
        z_acum = 0
        
        for j in range(individuals):
            string_bits = ''
                        
            if i == 0:
                val = False
                phenotype_a = 0.0
                phenotype_b = 0.0

                while not(val):
                    # Generación de la cadena de bit de acuerdo con el número de cromosomas
                    string_bits = generateBits(num_cromosmas)
                    # print(f'STRING: {string_bits}')
                    val, phenotype_a, phenotype_b = generatePhenotypesLineal( string_bits, linear, mja, mjb)

                # print('\n')
                z = computeZLinear(x, y, phenotype_a, phenotype_b)
                # Ingresamos los valores en nuestro diccionario
                linear_vectors[f'V{str(j+1)}'] = [ string_bits, phenotype_a, phenotype_b, z ]

            # Si no existe el vector V{index}
            else:
                if not( f'V{j+1}' in linear_vectors.keys()):
                    val = False
                    phenotype_a = 0.0
                    phenotype_b = 0.0
    
                    while not(val):
                        # Generación de la cadena de bit de acuerdo con el número de cromosomas
                        string_bits = generateVector(linear_vectors)
                        # print(f'STRING: {string_bits}')
                        val, phenotype_a, phenotype_b = generatePhenotypesLineal( string_bits, linear, mja, mjb)
                    
                    z = computeZLinear(x, y, phenotype_a, phenotype_b)
                    # Ingresamos los valores en nuestro diccionario
                    linear_vectors[f'V{str(j+1)}'] = [ string_bits, phenotype_a, phenotype_b, z ]

                else:
                    # Se obtiene el valor de Z de los vectores competentes seleccionados
                    z = linear_vectors[f'V{j+1}'][-1]                
                
            # Calculando la Z acumulada
            z_acum += z

        # print(f'LINEAR VECTORS: {linear_vectors}')

        # Proceso de selección para los vectores de la siguiente población
        aux = 0
        aux_2 = 0
        z_acums = []
        for _, vector in linear_vectors.items():
            aux = vector[-1] / z_acum
            aux_2 += (1-aux) / (individuals-1) 
            z_acums.append(aux_2)

        # Se geran números aletorios para cada individuo de la población
        random_numbers = [ random.random() for i in range(individuals) ]

        # print(f'Z ACUMULADAS: {z_acums}')
        # print(f'RANDOM NUMBERS: {random_numbers}')

        for rnd in random_numbers:
            for index, z_acum in enumerate(z_acums):
                if rnd < z_acum:
                    competent_vectors.append(index+1)
                    break
                
        # Vectores competentes para la siguiente iteración
        # print(f'\nCOMPETENT VECTORS: {competent_vectors}\n')
        occurences = { str(number) : competent_vectors.count(number) for number in competent_vectors }

        sort_ocurrences = sorted(occurences.items(), key = lambda x: x[1], reverse=True)
        # print(f'OCURRENCE: {sort_ocurrences}') 

        # Calculando el mejor vector para esta población
        best_vector = getBestZ( sort_ocurrences, linear_vectors )
        # print(f'BEST VECTOR OF POPULATION[{i+1}] = {best_vector}\n')

        competent_vectors = list(set(competent_vectors))
        competent_vectors.sort()
        # print(f'\nCOMPETENT VECTORS SIN REPETIR: {competent_vectors}\n')
        # Ingresamos los nuevos vectores competentes al diccionario
        aux_linear_vectors = dict()

        for index, vector in enumerate(competent_vectors):
            aux_linear_vectors[f'V{index+1}'] = linear_vectors[f'V{vector}']
        
        linear_vectors = aux_linear_vectors
        # print(f'\nNEXT POPULATION: {linear_vectors}\n')

    # Se retorna el mejor vector de cada población
    return best_vector
        
def computeCG( populations, individuals, linear, num_cromosmas, mja, mjb, mjc, x, y, flag ):

    linear_vectors = dict()
    best_vector = []
    
    for i in range(populations):
        competent_vectors = []
        z_acum = 0
        
        for j in range(individuals):
            string_bits = ''
                        
            if i == 0:
                val = False
                phenotype_a = 0.0
                phenotype_b = 0.0
                phenotype_c = 0.0

                while not(val):
                    # Generación de la cadena de bit de acuerdo con el número de cromosomas
                    string_bits = generateBits(num_cromosmas)
                    # print(f'STRING: {string_bits}')
                    val, phenotype_a, phenotype_b, phenotype_c = generatePhenotypesCG( string_bits, linear, mja, mjb, mjc)

                # print('\n')
                z = computeZCG(x, y, phenotype_a, phenotype_b, phenotype_c, flag)
                # Ingresamos los valores en nuestro diccionario
                linear_vectors[f'V{str(j+1)}'] = [ string_bits, phenotype_a, phenotype_b, phenotype_c, z ]

            # Si no existe el vector V{index}
            else:
                if not( f'V{j+1}' in linear_vectors.keys()):
                    val = False
                    phenotype_a = 0.0
                    phenotype_b = 0.0
                    phenotype_c = 0.0
    
                    while not(val):
                        # Generación de la cadena de bit de acuerdo con el número de cromosomas
                        string_bits = generateVector(linear_vectors)
                        # print(f'STRING: {string_bits}')
                        val, phenotype_a, phenotype_b, phenotype_c = generatePhenotypesCG( string_bits, linear, mja, mjb, mjc)
                    
                    z = computeZCG(x, y, phenotype_a, phenotype_b, phenotype_c, flag)
                    # Ingresamos los valores en nuestro diccionario
                    linear_vectors[f'V{str(j+1)}'] = [ string_bits, phenotype_a, phenotype_b, phenotype_c, z ]

                else:
                    # Se obtiene el valor de Z de los vectores competentes seleccionados
                    z = linear_vectors[f'V{j+1}'][-1]                
                
            # Calculando la Z acumulada
            z_acum += z

        # print(f'LINEAR VECTORS: {linear_vectors}')

        # Proceso de selección para los vectores de la siguiente población
        aux = 0
        aux_2 = 0
        z_acums = []
        for _, vector in linear_vectors.items():
            aux = vector[-1] / z_acum
            aux_2 += (1-aux) / (individuals-1) 
            z_acums.append(aux_2)

        # Se geran números aletorios para cada individuo de la población
        random_numbers = [ random.random() for i in range(individuals) ]

        # print(f'Z ACUMULADAS: {z_acums}')
        # print(f'RANDOM NUMBERS: {random_numbers}')

        for rnd in random_numbers:
            for index, z_acum in enumerate(z_acums):
                if rnd < z_acum:
                    competent_vectors.append(index+1)
                    break
                
        # Vectores competentes para la siguiente iteración
        # print(f'\nCOMPETENT VECTORS: {competent_vectors}\n')
        occurences = { str(number) : competent_vectors.count(number) for number in competent_vectors }

        sort_ocurrences = sorted(occurences.items(), key = lambda x: x[1], reverse=True)
        # print(f'OCURRENCE: {sort_ocurrences}') 

        # Calculando el mejor vector para esta población
        best_vector = getBestZ( sort_ocurrences, linear_vectors )
        # print(f'BEST VECTOR OF POPULATION[{i+1}] = {best_vector}\n')

        competent_vectors = list(set(competent_vectors))
        competent_vectors.sort()
        # print(f'\nCOMPETENT VECTORS SIN REPETIR: {competent_vectors}\n')
        # Ingresamos los nuevos vectores competentes al diccionario
        aux_linear_vectors = dict()

        for index, vector in enumerate(competent_vectors):
            aux_linear_vectors[f'V{index+1}'] = linear_vectors[f'V{vector}']
        
        linear_vectors = aux_linear_vectors
        # print(f'\nNEXT POPULATION: {linear_vectors}\n')

    # Se retorna el mejor vector de cada población
    return best_vector


"""" Esta función generá los vectores para la siguiente población haciendo mutaciones y cruces """
def generateVector( linear_vectors ):
    new_vector = ''
    # Seleccionamos un vector apto de manera aleatoria
    competent_vectors = [ vector[0] for vector in linear_vectors.values()]
    random.shuffle(competent_vectors)
    
    # Mutación y Cruce
    if len(competent_vectors) > 1:     
       method = 'mutation' if bool(random.getrandbits(1)) else 'crossover'
       if method == 'mutation':
           vector = competent_vectors[0]
        #    print(f'SELECTED VECTOR: {vector}')
           new_vector = mutation(vector)
        #    print(new_vector + '\n')
       else:
            vector_1 = competent_vectors[0]
            # print(f'SELECTED VECTOR 1: {vector_1}')
            vector_2 = competent_vectors[1]
            # print(f'SELECTED VECTOR 2: {vector_2}')
            new_vector = crossover(vector_1, vector_2)
            # print(new_vector + '\n') 
    # Mutación
    else:
        vector = competent_vectors[0]
        # print(f'SELECTED VECTOR: {vector}')
        new_vector = mutation(vector)
        # print(new_vector + '\n')

    return new_vector 

def mutation(bits):
    # print('MUTATION')
    # print(f'VECTOR: {bits}')
    bit = random.randint(0, len(bits)-1)
    # print(f'MUTATED BIT: {bit}')
    list_bits = list(bits)
    list_bits[bit] = '0' if list_bits[bit] == '1' else '1'
    return "".join(list_bits)

def crossover( bits_1, bits_2 ):
    # print('CROSSOVER')
    # print(f'VECTOR 1: {bits_1}')
    # print(f'VECTOR 2: {bits_2}')
    list_bits_1 = list(bits_1)
    list_bits_2 = list(bits_2)

    init_bits_1 = random.randint(0, len(bits_1)-2) # 4
    # print(f'INIT: {init_bits_1}') 
    shift_bits_1 = (len(bits_1)-1) - init_bits_1   # 3
    # print(f'SHIFT: {shift_bits_1}') 
    fin_bits_1 = random.randint(1, shift_bits_1)   # 2
    # print(f'FIN: {fin_bits_1}') 
    

    list_bits_2[init_bits_1:init_bits_1+fin_bits_1+1] = list_bits_1[init_bits_1:init_bits_1+fin_bits_1+1]
    return "".join(list_bits_2)

def computeZLinear(x, y, phenotype_a, phenotype_b ):
    # Z es el valor de la función evaluada en los puntos
    z = 0
    for k in range(len(x)):
        z += (y[k] - (x[k]*phenotype_a+phenotype_b)) ** 2    
    return z

def computeZCG(x, y, phenotype_a, phenotype_b, phenotype_c, flag ):
    # Z es el valor de la función evaluada en los puntos
    z = 0
    if flag:
        for k in range(len(x)):
            z += (y[k] - ( (phenotype_a * x[k] ** 2 ) + (phenotype_b * x[k]) + phenotype_c )) ** 2
    else:
        for k in range(len(x)):
            z += (y[k] - ( phenotype_a * math.e ** ( -phenotype_b * (x[k] - phenotype_c)**2 )))** 2
    return z


def getBestZ( ocurrences, linear_vectors ):
    # Seleccionamos los mejores vectores de cada población
    best_vectors = []
    # Obtenemos el valor más grande
    threshold = ocurrences[0][1]

    for key, value in ocurrences:
        if value == threshold:
            best_vectors.append( linear_vectors[f'V{key}'][1:] ) 
    
    if len(best_vectors) > 1:
        best_vector = best_vectors[0]
        minimum = best_vectors[0][-1] 
        # Seleccionamos el menor valor
        for vector in best_vectors:
            if  vector[-1] < minimum:
                minimum = vector[-1]
                best_vector = vector
        
        best_vector[-1] = -best_vector[-1]
        return best_vector

    else:
        best_vectors[0][-1] = - best_vectors[0][-1]
        return best_vectors[0]    

