    #!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import subprocess
import argparse
import logging
import fnmatch
import pandas as pd
from os.path import join # Edit 

#para usar este script utilizando otros genes hay que cambiar en la linea 39 directorio_script+".."
#ejemplo (gen_p97.txt es mi archivo de genes, debe contener separado por columnas: las coordenadas, la orientacion y el codigo Rv)
#como en este caso: 522347  524533  R  Rv0435c_Rv0435c
#        subprocess.call("""perl """+ directorio_script+ """/findGeneDeletions2.pl """+gnumber+".coverage.all.pos "+ directorio_script+"/gen_p97.txt 10 > "+gnumber+".deletions", shell = True)


def read_fileGnumbers(fichero):
    lista_gnumbers = list()
    with open(fichero,"r") as g_file:
        for gnumber in g_file:
            if "G" in gnumber:
                gnumber = gnumber.replace("\n","")
                lista_gnumbers.append(gnumber)
    logging.warning(lista_gnumbers)
    return lista_gnumbers
def generate_deletions_file(gnumbers_list,ruta_original,directorio_script):
    for gnumber in gnumbers_list:
        os.chdir(ruta_original + '/' + gnumber[0:3] + '/' + gnumber[3:5] + '/' + gnumber[5] + '/')
        subprocess.call("""gunzip """+gnumber+".all.pos.vcf.gz", shell = True)
        subprocess.call("""awk -F '[\t;=]' '{if ($0 !~/#/) {print $2"\t"$9}}' """+gnumber+".all.pos.vcf > "+gnumber+".coverage", shell = True)
        subprocess.call("""perl """+ directorio_script+ """/built_total_coverage.pl """+gnumber+".coverage > "+gnumber+".total.coverage", shell = True)

        listOfFiles = os.listdir(ruta_original + '/' + gnumber[0:3] + '/' + gnumber[3:5] + '/' + gnumber[5] + '/')
        ''' cambiado chan pq algunos ya tenían el coverage calculado por paula
        pattern = "*.coverage.all.pos"
        for entry in listOfFiles:
            if fnmatch.fnmatch(entry, pattern):
                logging.warning("The coverage.all.pos file already exists, i won't overwrite it.")
                continue
            else:
                subprocess.call("""awk '{print $2}'  """+gnumber+".total.coverage > "+gnumber+".coverage.all.pos", shell = True)
        '''
        subprocess.call("""awk '{print $2}'  """+gnumber+".total.coverage > "+gnumber+".coverage.all.pos", shell = True)
        subprocess.call("""rm """+gnumber+".coverage", shell = True)
        subprocess.call("""rm """+gnumber+".total.coverage", shell = True)
        subprocess.call("""gzip """+gnumber+".all.pos.vcf", shell = True)
        subprocess.call("""perl """+ directorio_script+ """/findGeneDeletions2.pl """+gnumber+".coverage.all.pos "+ directorio_script+"/allGenesMTBC.txt 10 > "+gnumber+".deletions", shell = True)
        

def main():
    mapsPath = '/storage/PGO/data/mtb/mappings/v1/' ### Edit
    resultsPath = '/storage/PGO/results/mtb/deletions_allGenes/' ###Edit
    directorio_script = os.getcwd()
    parser = argparse.ArgumentParser(description = 'Download fastq given a B number and create requiered directory.')
    parser.add_argument('-f', dest = 'file', help = 'File with numbers.', required = True)
    parser.add_argument('-d', dest = 'directory', help = 'Directory where B numbers are stored.', required = True)

    args = parser.parse_args()
    ruta_original = args.directory
    gnumbers_list = read_fileGnumbers(args.file)
    generate_deletions_file(gnumbers_list,ruta_original,directorio_script)
    os.chdir(directorio_script)

    #cambiado chan pq no existe el directorio
    #human_file_T = open("genes_"+args.file,"w")
    human_file_T = open(args.file,"w")
    deletionsDF = pd.DataFrame()# Edit
    for gnumber in gnumbers_list:
        line_human_T = gnumber+"\t"

        os.chdir(ruta_original + '/' + gnumber[0:3] + '/' + gnumber[3:5] + '/' + gnumber[5] + '/')
        nombre = gnumber+".deletions"
        percList = []# Edit
        geneList = []# Edit
        with open(gnumber+".deletions","r") as g_file:
            for line in g_file:
                line = line.replace("\n","")
                if "Deletion:" in line:
                    line_split = line.split(" ")
                    linea = line_split[1] + " " + line_split[2]
                    line_human_T+= linea +"\t"
                    perc = float(line_split[2])# Edit
                    percList.append(perc)# Edit
                    geneList.append(line_split[1])# Edit
        percList = pd.Series(percList)# Edit
        deletionsDF = pd.concat([deletionsDF, percList], axis = 1)# Edit       
        #subprocess.call("""rm """+gnumber+".deletions", shell = True)
        #logging.warning("the .deletions file was removed")
        os.chdir(directorio_script)
        line_human_T += "\n"

        human_file_T.write(line_human_T)
        
        gNumTxt = pd.concat([pd.Series(geneList), percList], axis = 1)### Edit 
        pathOut = os.path.join(mapsPath, gnumber[0:3], gnumber[3:5], gnumber[5]) ### Edit
        pathOut = pathOut + '/' + gnumber + '_dels.txt' ### Edit
        gNumTxt.to_csv(pathOut, header=None, index=None, sep='\t', mode='a') ### Edit
        print(gnumber + " deletions saved at " + pathOut + ".") ### Edit
        
    deletionsDF.columns = gnumbers_list# Edit
    deletionsDF.index = geneList# Edit
    endName = str(args.file)# Edit
    endName = endName.split('/')[6]# Edit
    endName = endName.split('.')[0]# Edit
    deletionsDF.to_csv(resultsPath + 'deletionsDF_' + endName + '.csv')# Edit   
    print("deletionsDF_" + endName + ".csv saved at " + resultsPath + ".") ### Edit
    human_file_T.close()
    

if __name__ == "__main__":
    #print(__doc__) #imprimimos la documentacion
    main()