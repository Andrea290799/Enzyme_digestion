#!/usr/bin/env python
# -*- coding: latin1 -*-
# Joaquin Giner Lamia 2020

#importamos modulos
import sys
import re


#----- functions -------


# Function for error controlling 

def help():
	
	message = (

"""#_Enzyme digestion script_# 

This program perform a double DNA digestion using the 
restriction enzymes abaI and abaII. The program obtain the number of DNA fragments and
their length after this double enzyme digestion


Usage: python enzyme_digestion.py [-h] <fasta_file>

	
	-h: help
	fasta_file: a fasta format file containing the DNA sequence for digestion

	"""
			)

	print(message)
	sys.exit()


def warning():
	
	message = (

"""Error: few parameters
Usage: python enzyme_digestion.py [-h] <fasta_file>

	
	-h: help
	fasta_file: a fasta format file containing the DNA sequence for digestion

	"""

			)
	
	print(message)
	sys.exit()


def open_input_file(input_file):
	
	"""
	abre el archivo de entrada, lo lee y elimina los saltos de linea
	"""
	
	dna = open(input_file).read().rstrip("\n")

	return dna


def is_fasta(dna):
	"""
	controla que el archivo de entrada sea fasta en caso negativo aborta el
	script
	"""

	result = False

	for line in dna:
		if line.startswith(">"):
			result = True
			break
	if result == False:
		sys.exit("\nError:provided DNA sequence file is not in fasta format!")
		


def input_control():
	"""
	Controla todas las excecipiones que puede encontrar en los paramtros de entrada
	que no haya parametros, que este el paramtro -h, que no se pueda abrir el archivo
	fasta o que no sea un archivo fasta. Todo esto controlado por otras funciones que en
	caso de no cumplir las condiciones exigidas provocaran la salida del script con la
	funcion sys.exit()
	"""

	arguments = sys.argv
	
	# si el numero de argumentos es menor de 2, es decir no hay argumentos mas alla del propio script
	# provocamos un warning
	if len(arguments) < 2:
		warning()

	else:
		# iteramos sobre los argumentos
		for argument in arguments:
			
			if argument == "-h":
				help()

			else:	
				# si el argumento no es -h intentamos leerlo
				try:
					dna_file = open_input_file(argument)	

				# si da fallo la lectura soltamos un error y terminamos ejecucion
				except:
					print("\nError: There is a problem with your fasta file ")
					sys.exit()

	# si todo va bien hasta aqui comprobamos si es fasta				
	is_fasta(dna_file)				
					
	return dna_file




# Function for digestion computing


def enzymes_cuts(dna_file):

	"""
	genera una lista que contiene todos las posiciones de los cortes con las enzimas 
	de restricción.
	"""
	
	# inicialiazmos la variable que contiene la lista de cortes añadiendo el punto cero,
	# esto es porque la secuencia de ADN es lineal y no circular
	all_cuts = [0]


	# cortamos con AbcI, sumando la posicion exacta del corte en este caso +3
	for match in re.finditer(r"A[ATGC]TAAT", dna_file):
		all_cuts.append(match.start() + 3)


	# cotamos con AbcII, sumando la posicion exacta del corte en este caso +4	
	for match in re.finditer(r"GC[AG][AT]TG", dna_file):
		all_cuts.append(match.start() + 4)


	# añadimos a la lista el final de la secuencia de ADN.	
	all_cuts.append(len(dna_file))

	# ordenamos los sitios de cortes a lo largo de la secuencia de menor a mayor
	sorted_cuts = sorted(all_cuts)
	
	# imprimimos todos los cortes
	print("\nRestriction enzymes cut in:"+str(sorted_cuts[1:len(sorted_cuts)-1]))

	return sorted_cuts



def compute_fragment_and_print_result(sorted_cuts):
	
	"""
	itera sobre la lista de sorted_cuts y va restando al valor x x-1, para obtener 
	la longitud de los fragmentos de ADN. El numero total de fragmentos lo calula
	contando el numero de iteraciones.
	La función muestra el numero total de fragmentos en pantalla y genera la lista 
	de fragmentos que sera utilizada por la funcion graphic

	"""

	# inicializamos variable que acumula el numero de fragmentos
	fragments = 0

	# inicializamos lista donde almacenaremos fragmentos
	fragments_list = []

	# iteramos sobre la lista sorted_cuts para calcular los fragmentos de ADN.
	for i in range(1,len(sorted_cuts)):
		cut_position = sorted_cuts[i]
		previous_cut_position = sorted_cuts[i - 1]
		fragment_size = cut_position - previous_cut_position
		#print("fragment size is: "+ str(fragment_size))

		fragments_list.append(fragment_size)
		fragments +=1

	print("Total DNA fragments after digestion: "+ str(fragments))


	return fragments_list


# Function for graphication

def graphic_restriction(fragments_list):
	
	"""
	Esta funcion grafica los fragmentos cortados de DNA utilizando
	la lista fragment list
	"""

	# DNA sizes:
	small = "##    "
	medium = "####  "
	big = "######"


	# ordenamos la lista de fragmentos
	length_fragment = sorted(fragments_list, reverse=True)
	
	print(("\nFragments prediction for: {}\n").format(sys.argv[1]))	

	# iteramos sobre la lista de fragmentos de ADN
	for i in range(0,len(length_fragment)):
		
		# asignamos un tamaño de banda de restricción en funcion del tamaño 
		dna_length = length_fragment[i]
		if dna_length < 200:
			shape = small
		elif 200 < dna_length < 500:
			shape = medium 	
		else:
			shape = big

		print(("\n   {}   {}bp").format(shape, str(dna_length)))

	print("\n")






def main():

	#obtenemos el fasta input y controlamos errores:
	dna_file = input_control()

	#si todo ha ido bien digerimos
	sorted_cuts = enzymes_cuts(dna_file)

	#calculamos fragmentos de ADN y su longitud 
	fragments_list = compute_fragment_and_print_result(sorted_cuts)
		
	#imprimimos los fragmentos en pantalla
	graphic_restriction(fragments_list)



# Solo se ejecutará la funcion main si el codigo es ejecutado como script, si es llamado
# como modulo por otro script la funcion main no se ejecutará. Es decir si sys.argv[0] 
# es igual al nombre de este script ejecuta la función main(input_file). 
if __name__ == '__main__':
	main()









