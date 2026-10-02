#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Noemí Toro Barrios | TFM - Revisión de herramientas y recursos para la anotación y priorización de variantes en regiones no codificantes

#Funciones auxiliares para parsear ClinVar, asignar el número de estrellas,
#clasificar las variantes y aplicar los filtros de evidencia y consecuencia.


def parse_info(info_str):
    #Convierto la información en un diccionario para poder llamar a cada campo de forma sencilla. 
    
    
    if info_str == ".": #Si el campo está vacío, no me devuelvas nada
        return {}

    info_dict = {} #Creo el diccionario vacío

    #Separo los campos que aparecen en el VCF por ; 

    for part in info_str.split(";"): 

        if '=' in part: 

            k, v = part.split("=",1) 
            info_dict[k] = v #Separa la anotación en el diccionario para que cuando llamemos a MC nos de la consecuencia y así con todos los campos

            #En caso de que haya alguna anotación que no esté separado por igual, añado esa línea como TRUE       

        else:
            
            info_dict[part] = "NO HAY VALOR"

    return info_dict


def calcular_estrellas(revstat):

    #Filtrar por 2 o más entrellas. 

    if not revstat or revstat==".": return 0 #Si está vacío, me da 0 estrellas directamente. No quiero esa variante. 

    #Filtrado de estrella según: https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/

    if "practice_guideline" in revstat: return 4  

    if "reviewed_by_expert_panel" in revstat: return 3

    if "criteria_provided" in revstat and "multiple_submitters" in revstat and "no_conflicts" in revstat: return 2

    if "criteria_provided" in revstat and "single_submitter" in revstat: return 1

    return 0 #Si el campo no cumple con nada de lo anterior, no nos fiamos. 



def clasificar_variante(clnsig, clnsig_conf):

    #Clasifico las variantes en 5 grupos: benigno/probablemente benigno, VUS, conflicting, patogénico/probablemente patogénico y otras

    #Normalizo todo el campo a minúscula para no tener problemas a la hora de llamar una categoría

    if clnsig:
        sig= clnsig.lower()
    else:
        sig= ""

    if clnsig_conf:
        conf= clnsig_conf.lower()
    else:
        conf= ""

    #Automáticamente añado la variante a cateogria otras

    categoria = "Otras"

    #Como luego me quedaré solo con aquellas conflicting que tengan una entrada patogénica, reviso si es patogénica o no ahora:

    patogenico = False

    if "conflicting" in sig: #Variante conflicting:

        categoria = "Conflicting"

        if "pathogenic" in conf:

            patogenico = True

    elif "pathogenic" in sig:

        categoria = "Patogénica/Probablemente patogénica"
    elif "benign" in sig:

        categoria = "PBenigna/Probablemente benigna"

    elif "uncertain" in sig:

        categoria = "VUS"

    return categoria, patogenico

def filtrar(estrellas, patogenico):

    #Filtrar por número de estrellas o conflicting con entradas patogénica

    if estrellas >= 2:
        return True

    if patogenico:

        return True

    return False



    #Filtrar solo por consecuencia no codificante. Si tiene una consecuencia codificante descartamos. 


PALABRAS_CODIFICANTES = ["missense", "nonsense", "synonymous", "frameshift", "inframe", "stop_gained", "stop_lost", "start_lost", "protein_altering", "coding_sequence"]

def no_codificante(mc_raw):
    if not mc_raw or mc_raw == "." or mc_raw == "NO HAY VALOR":
        return False

    consecuencias = mc_raw.split(".")

    #Elimino las codificantes. 

    for c in consecuencias:
        c_minuscula = c.lower() #ponemos todo en minusculas para que no tengamos errores al comparar las palabras

        for prohibida in PALABRAS_CODIFICANTES:
            if prohibida in c_minuscula:

                return False
            
    return True

 
