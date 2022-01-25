# Instale Bioconductor version 3.14 (dejo como por si acaso):
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

## Primero cargue las librerias:
library(BiocManager)
library(Biostrings)

### Secuencia correspondiente de aminoácidos:
sc <- readRNAStringSet ("secuencias.fasta")
aa <- translate (sc)
aa

##############################
# Plataforma Rosalind:

#### Ejercicio 1: Complementing a Strand of DNA ####
### Enlace: https://rosalind.info/problems/revc/

## Con librerias especializadas: Biostrings
adn <- DNAString ("AAAACCCGGT")
reverseComplement (adn)
# ejemplo en Rosalind

gen_casei <- readDNAStringSet ("gen_shirota.txt")
reverseComplement (gen_casei)
# con una secuencia de mas de 1000 bp
# en este caso: el gen mntH3 de Lactobacillus casei

# Use la libreria Biostrings, la cual cargue arriba
# Primero con la funcion DNAString y/p readDNAStringSet especifique la cadena que tome del problema de Rosalind... 
# ...es una cadena de ADN y lo converti ademas en un objeto llamado adn y gen_casei
# Despues con la funcion reverseComplement realice todo, porque tal funcion ya realiza lo que queriamos...
# ...lo cual ver la hebra antisentido de ADN, de derecha a izquierda.

## Sin libreria especializada:
sec_adn <- c ("A","A","A","A","C","C","C","G","G","T")

rep_1 <- replace (sec_adn, sec_adn == "A", "t")
rep_2 <- replace (rep_1, rep_1 == "T", "a")
rep_3 <- replace (rep_2, rep_2 == "C", "g")
rep_4 <- replace (rep_3, rep_3 == "G", "c")
rep_4 #Para observar el resultado de el cambio de todas las bases
rever<- rev(rep_4)
rever # Resultado

# Primero coloque mi secuencia en un objeto llamado sec_adn, en donde separe su composicion en comillas...
#... esto es con el proposito de que R lea cada letra por separado.
# Despues con la funcion replace reemplace A-T y C-G, sin embargo, como R puede reconocer entre...
# ... minusculas y mayusculas, decidi cambiarla por la base pero en minuscula (es decir, A-t).
# Esto lo realice porque si lo hacia con mayuscula en rep_1 cambiaria todas las A?s por T?s y debido...
# ... a que esta encadenado/escalonado cuando llegara a rep_2 regresaria las T?s que ya habia cambiado a A?s. Y eso...
#... tambien pasaria con G y C.
# Asi que despues de haber cambiado todas las letras, para poder realizar el reverso use la funcion rev...
# ... que tal cual esa es su funcion.
# Nota: Tambien encontre que en lugar de usar la funcion replace se puede usar chartr.

##### Ejercicio 2: Counting DNA Nucleotides ####
### Enlace: https://rosalind.info/problems/dna/

## Con librerias especializadas: Biostrings

nuc <- DNAString ("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC")
oligonucleotideFrequency(nuc, 1)
# ejemplo en Rosalind

gen_casei <- readDNAStringSet ("gen_shirota.txt")
oligonucleotideFrequency(gen_casei, 1)
# con una secuencia de mas de 1000 bp
# en este caso: el gen mntH3 de Lactobacillus casei

# Use la libreria Biostrings, la cual cargue arriba
# Primero con la funcion DNAString y/p readDNAStringSet especifique la cadena que tome del problema de Rosalind... 
# ...es una cadena de ADN y lo converti ademas en un objeto llamado gen_casei y nuc.
# Despues con la funcion oligonucleotideFrequency tal cual cuenta los nucleotidos, por lo que no debemos agregar mas.
# El 1 de al final especifica las combionaciones de las bases, en este caso las queremos contgar por separado, por eso 1.

## Sin libreria especializada:

nchar (gen_casei)

nucl_a <- grepRaw ("A",gen_casei, all=T)
total_ade <- length (nucl_a)
total_ade # Primer resultado: adenina

nucl_c <- grepRaw ("C",gen_casei, all=T)
total_cit <- length (nucl_c)
total_cit # Segundo resultado: citocina

nucl_g <- grepRaw ("G",gen_casei, all=T) 
total_gua <- length (nucl_g)
total_gua # Tercer resultado: guanina

nucl_t <- grepRaw ("T",gen_casei, all=T)
total_tim <- length (nucl_t)
total_tim # Cuarto resultado: timina

# Primero con la funcion nchar podemos observar el numero de caracteres que hay.
# Con la funcion grepRaw busca la base que se le esta pidiendo, en estos casos es la base que se encuentra...
#... en los parentesis y que esta entre comillas (ej. "A").
# Despues con length cuenta lo que se encontro antes con grep y nos muestra ahora si el numero de bases que encontro.
# Tambien all es para que lo lea en todo, no solo en uno.
