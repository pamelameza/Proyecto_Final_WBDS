# Proyecto_Final_WBDS

# **Women in Bioinformatics and Data Science**

## <font color=purple>**2023 1st WBDS LA Camp "Introduction to Bioinformatics and Data Science"** 
Proyecto Final

Elaborado por Pamela Denisse Meza Cruz

Objetivo: Búsqueda de secuencias codificantes de proteínas asociadas a la resistencia a los antibióticos
  
# **0. Preparación del entorno**

```python
!pip3 install matplotlib
!pip3 install pandas
!pip3 install pycirclize
!pip3 install pyrodigal
!pip3 install requests
!pip3 install seaborn
!pip3 install biopython
```
```python
import matplotlib.pyplot as plt
import numpy  as np
import pandas as pd
import sys
import pyrodigal
import requests
import yaml
import seaborn as sns
import subprocess
import re
!pip3 install --upgrade Biopython
!pip3 install pycirclize
from Bio import Entrez
from Bio import SeqIO
Entrez.email   = "pamelamezac13@gmail.com"
Entrez.api_key = "e98f23ed08456bc278dcdbc65c02a160a708"
from pathlib import Path
from time import sleep
from typing import TextIO, Union
from io                 import StringIO
from matplotlib.patches import Patch
from pycirclize         import Circos
from pycirclize.parser  import Gff
from requests.adapters  import HTTPAdapter, Retry
```
  
# **1. Obtención de una secuencia genómica**
  
### Búsqueda en el sitio web de NCBI del organismo de interés, en este proyecto se seleccionó *Klebsiella*, un género de bacterias considerado como uno de los principales patógenos que generan resistencia a los antibióticos

- #### Para este proyecto se utilizó el ensamble [ASM763225v1](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_007632255.1/), y se buscaron genes que están implicados en la resistencia a los antibióticos en el cromosoma con número de acceso [CP041925.1](https://www.ncbi.nlm.nih.gov/nuccore/CP041925.1/)

### Descargar el genoma de nuestro interés en formato genbank.
```python
accession = "CP041925.1"
genome = Entrez.efetch(db="nucleotide",
                       id=accession,
                       format="gb",
                       rettype="text")
record = SeqIO.read(genome, "genbank")
genome_length = len(record.seq)
  
genome_length
5183422
```  
# **2. Predicción de genes usando pyrodigal**
  
### Encontrar los genes codificantes de proteínas en el genoma procariótico
```python
orf_finder = pyrodigal.OrfFinder()
orf_finder.train(bytes(record.seq))
orf_genes  = orf_finder.find_genes(bytes(record.seq))
```
### Almacenar las secuencias aminoacídicas de los genes predichos en un nuevo archivo `CP041925.1.faa`

> ### Es importante especificar un prefijo para identificar nuestras secuencias aminoacídicas. En este caso estamos usando el genoma de *Klebsiella*, por lo que le pondremos Kleb
```python  
aa_file = accession + ".faa"
prefix  = "Kleb"
with open(aa_file, "w") as orf_gene:
    orf_genes.write_translations(orf_gene,sequence_id=prefix)
 ``` 
### Almacenar las coordenadas de los genes predichos en un nuevo archivo `CP041925.1.gff`
```python  
gff_file = accession + ".gff"
prefix  = "Kleb"
with open(gff_file, "w") as orf_gene:
    orf_genes.write_gff(orf_gene,sequence_id=prefix)
```
! head CP041925.1.gff
  ![Captura de pantalla 2023-03-16 a la(s) 1 31 24 p m](https://user-images.githubusercontent.com/117956356/225732878-47aa0bad-35b7-4496-9041-d3c077205e4e.png)

- ### El archivo `CP041925.1.faa` será util para comparar un set de secuencias curadas contra las predicciones de pyrodigal
- ### El archivo `CP041925.1.gff` será util para visualizar las predicciones realizadas
  
  
# **3. Obtención de un set de secuencias de referencia**

## Utilizar la API de UniProt

### Descargar secuencias de UniProt en el objeto `uniprot_ref_seqs`
```python  
uniprot_api_url  = "https://rest.uniprot.org/uniprotkb/stream"
uniprot_api_args = {"compressed" : "false",
                    "format"     : "fasta",
                    "query"      : "(antibiotic resistance) AND (reviewed:true)"}
uniprot_ref_seqs = requests.get(uniprot_api_url,params=uniprot_api_args).text
 ``` 
### Pasar las secuencias en formato fasta a un archivo (`uniprot_sequences.fasta`) en nuestro disco duro
``` python  
uniprot_seqs_file = open("uniprot_sequences.fasta", "wt")
uniprot_seqs_file.write(uniprot_ref_seqs)
uniprot_seqs_file.close()
```  
! head uniprot_sequences.fasta
![Captura de pantalla 2023-03-16 a la(s) 1 32 57 p m](https://user-images.githubusercontent.com/117956356/225733163-32abccc6-48e8-4f4c-9859-ac2cb6462fc1.png)

  
# **4. Creación de una base de datos tipo BLAST**
  
## A partir de las secuencias aminoacídicas de las predicciones de pyrodigal, crear una base de datos tipo BLAST que se usará para recopilar secuencias similares a las secuencias obtenidas en UniProt
``` python
makeblastdb_path = "/Users/Pamela/Downloads/ncbi-blast-2.7.1+/bin/makeblastdb"
makeblastdb_command = [makeblastdb_path,'-in',aa_file,'-dbtype','prot']
subprocess.call(makeblastdb_command)
 ``` 
![Captura de pantalla 2023-03-16 a la(s) 1 33 55 p m](https://user-images.githubusercontent.com/117956356/225733310-f93b0011-f430-4b75-916f-a888de275ac7.png)

# **5. Obtención de secuencias de interés en el genoma analizado**
  
## Obtención de dos sets de proteínas:

- ### 4760 Proteínas en el genoma
- ### 2714 Proteínas de referencia

## Comparar ambos sets para ver que secuencias en nuestro genoma podrían estar involucradas en la resistencia a antibioticos  
  
## Llamar a BLAST y obtener una tabla con los resultados de la comparación

### Comparar las secuencias de uniprot contra las secuencias del genoma  
```python  
blastp_path       = "/Users/Pamela/Downloads/ncbi-blast-2.7.1+/bin/blastp"
blastp_out_format = "6 qseqid sseqid qlen slen qstart sstart qend send score evalue length positive"
blastp_out_file   = accession + ".blast.tsv"
blastp_command    = [blastp_path,
                     "-db",          aa_file,
                     "-query",       "uniprot_sequences.fasta",
                     "-evalue",      "1e-9",
                     "-out",         blastp_out_file,
                     "-outfmt",      blastp_out_format,
                     "-num_threads", "4"]
subprocess.call(blastp_command)
```  
! head CP041925.1.blast.tsv
![Captura de pantalla 2023-03-16 a la(s) 1 35 18 p m](https://user-images.githubusercontent.com/117956356/225733533-9c3e0837-98fe-43b2-97fc-88503bdfcfe5.png)

## Comparar las secuencias obtenidas del genoma contra las secuencias de uniprot 
```python  
makeblastdb_path = "/Users/Pamela/Downloads/ncbi-blast-2.7.1+/bin/makeblastdb"
makeblastdb_command = [makeblastdb_path,'-in',"uniprot_sequences.fasta",'-dbtype','prot']
subprocess.call(makeblastdb_command)
blastp_path      = "/Users/Pamela/Downloads/ncbi-blast-2.7.1+/bin/blastp"
blast_out_format = "6 qseqid sseqid qlen slen qstart sstart qend send score evalue length positive"
blast_out_file   = "uniprot_sequences.blast.tsv"
blastp_command   = [blastp_path,
                    "-db",          "uniprot_sequences.fasta",
                    "-query",       aa_file,
                    "-evalue",      "1e-9",
                    "-out",         blast_out_file,
                    "-outfmt",      blast_out_format,
                    "-num_threads", "4"]
subprocess.call(blastp_command)
```
![Captura de pantalla 2023-03-16 a la(s) 1 36 15 p m](https://user-images.githubusercontent.com/117956356/225733748-ef7779b1-4107-4877-9d7c-163d536a2056.png)

# **6. Examinación de los resultados de la búsqueda tipo BLAST**
 
### Los resultados de la búsqueda tipo BLAST están en el archivo `CP041925.1.blast.tsv`, sin embargo dicha tabla no contiene los nombres de las columnas

### Podemos importar esta tabla a un dataframe de pandas y asignar los nombres de las columnas usando la variable `blastp_out_format` que definimos anteriormente
 ```python
blastp_column_names = blastp_out_format.split(" ")[1:]
blastp_df = pd.read_csv(blastp_out_file,sep="\t",names=blastp_column_names)
blastp_df
```  
![Captura de pantalla 2023-03-16 a la(s) 1 37 11 p m](https://user-images.githubusercontent.com/117956356/225733948-1db7616b-d983-4a28-b2c7-09363c7e45cc.png)

### En el genoma de interés encontramos 625 proteínas potencialmente asociadas con biosíntesis de antibioticos, no obstante, no todas las secuencias podrían ser de nuestro interés
  ```python
candidate_genes=blastp_df["sseqid"].unique().tolist()
len(candidate_genes)
625
  ```
## 6.1. Visualización preliminar de los datos
  
### Con la biblioteca pyCirclize podemos visualizar los genes que fueron  identificados en el paso anterior

### 6.1.1 Transformar nuestro archivo gff en un dataframe
```python  
gff_columns     = ["chr","source","feature_type","start","end","score","strand","phase","info"]
gff_df          = pd.read_csv(gff_file,sep="\t",comment="#",header=None,names=gff_columns)
gff_df["start"] = gff_df["start"].astype(int)
gff_df["end"]   = gff_df["end"].astype(int)
gff_df
  ```
![Captura de pantalla 2023-03-16 a la(s) 1 38 32 p m](https://user-images.githubusercontent.com/117956356/225734207-9109c38f-88ed-4145-9f40-8d12534c2341.png)

### 6.1.2 Obtención de información adicional del dataframe `gff_df`
```python
def get_gff_info(info_str):
    out_dict = {}
    info_arr = info_str.split(";")
    for line in info_arr:
        if "=" in line:
            line_arr    = line.split("=")
            field_name  = line_arr[0]
            field_value = line_arr[1]
            out_dict[field_name] = field_value
    return out_dict
```  
  ```python
gff_df["annotation"] = gff_df["info"].apply(lambda x: get_gff_info(x))
gff_df
```
  
### 6.1.3 Filtrado de datos para incluir solamente los genes identificados como asociados a la resistencia a antibióticos
```python  
gff_df["candidate"] = gff_df["annotation"].apply(lambda x: "include" if x["ID"] in candidate_genes else "exclude")
gff_df
 ``` 
  
### 6.1.4 Almacenar el resultado en un nuevo archivo gff para que pyCirclize visualice unicamente los genes de interés
 ```python 
candidate_df = gff_df.copy()
candidate_df = candidate_df[candidate_df["candidate"]=="include"][gff_columns]
candidate_df.to_csv("candidates.gff",sep="\t",header=False,index=False)
 ``` 
! head candidates.gff
![Captura de pantalla 2023-03-16 a la(s) 1 41 31 p m](https://user-images.githubusercontent.com/117956356/225734885-30a41afe-29b7-4a0a-9674-3470ffa1f647.png)

### 6.1.5 Visualización de los datos con pyCirclize

### Construir distintos objetos para obtener un mapa circular que nos permitirá identificar manualmente potenciales operones en el genoma de Klebsiella.
 ```python 
circos = Circos(sectors={accession: genome_length})
circos.text("Klebsiella aerogenes")
circos_gff = Gff(gff_file="candidates.gff")
sector = circos.get_sector(accession)
sector = circos.sectors[0]
cds_track = sector.add_track((80, 100))
cds_track.axis(fc="#EEEEEE", ec="none")
cds_track.genomic_features(circos_gff.extract_features("CDS", target_strand =  1), r_lim=(90, 100),fc="red" )
cds_track.genomic_features(circos_gff.extract_features("CDS", target_strand = -1), r_lim=(80,  90),fc="blue")
pos_list, labels = [], []
cds_track.xticks_by_interval(
    interval=500000,
    label_formatter=lambda label_value: f"{label_value/ 1000000:.1f} Mb",
    label_orientation="vertical")
fig = circos.plotfig().set_figwidth(5)
  ```
![Captura de pantalla 2023-03-16 a la(s) 1 42 16 p m](https://user-images.githubusercontent.com/117956356/225735079-a5929eff-dd0f-4b20-ae04-173a8f19b747.png)

### 6.1.6 Visualización de los datos con seaborn

### Si queremos realizar una examinación más exhaustiva, podemos usar una serie de swarmplots en donde podemos comparar las posiciones de los genes en el dataframe completo, separando por categorias ("genes candidatos" vs "genes no candidatos") y por cadena ("+" vs "-")

### Con esta aproximación podemos identificar un par de operones enriquecido en genes candidatos, siendo el operón ubicado en la cadena negativa entre 0.0 Mbp y 0.5 Mbp
```python  
num_bins = 25
counter_1 = 0
counter_2 = 0
fig, axes = plt.subplots(5,5,figsize=(30,30))
bin_len  = (genome_length - (genome_length % (num_bins - 1))) / (num_bins)
for bin_num in range(num_bins):
    start_pos = bin_num * bin_len
    end_pos   = (bin_num + 1) * bin_len
    mb_df = gff_df.copy()
    mb_df = mb_df[(mb_df["start"]>start_pos) & (mb_df["end"]<=end_pos)]
    sns.swarmplot(ax = axes[counter_1,counter_2],data = mb_df,y="candidate",x="start",hue="strand",dodge=True,order=["exclude","include"],hue_order=["+","-"])
    axes[counter_1,counter_2].set(ylabel=None)
    counter_2 += 1
    if (counter_2%5 == 0):
        counter_2 = 0
        counter_1 += 1
plt.show()
```                                                                      
![Captura de pantalla 2023-03-16 a la(s) 1 43 08 p m](https://user-images.githubusercontent.com/117956356/225735286-d0e0cf26-d61f-4ea4-8522-51403a108dd9.png)
                                                   
                                                                       
# **7. Examinación a detalle del operón seleccionado**
                                                                       
### Para obtener aún más información acerca del operón que seleccionamos, podemos analizar las secuencias aminoacídicas de los genes de dicha región a través del servicio de búsqueda de dominios conservados de InterProScan

### Lo primero que debemos hacer es obtener los IDs de los genes presentes en dicho operón                                                                 ```python      
operon_df = gff_df.copy()
operon_df = operon_df[(operon_df["start"]     >= 400000) &
                      (operon_df["end"]       <= 500000) &
                      (operon_df["strand"]    == "+")     &
                      (operon_df["candidate"] == "include")]
operon_df.reset_index(drop=True, inplace=True)
                                                       
len(operon_df)
33                                                       
```
```python                                                       
operon_gene_list = []
for index in operon_df.index.tolist():
    gene_id = operon_df["annotation"][index]["ID"]
    operon_gene_list.append(gene_id)                                                       
  
operon_gene_list
 ```                                                      
### Posteriormente, construir un *string* que contendrá las secuencias aminoacídicas de los genes de interés en formato fasta

### Enviar el *string* al servicio web de InterProScan para buscar dominios conservados en nuestras proteínas

> ### Antes de correr el servicio web, tenemos que eliminar los asteriscos de nuestras secuencias
```python  
query_str = ""
for record in SeqIO.parse(aa_file, "fasta"):
    seq_id  = record.id
    if(seq_id in operon_gene_list):
        seq_str = str(record.seq)
        query_str+=">"+seq_id+"\n"+seq_str+"\n"
query_str = query_str.replace("*","")
```  
### El proceso de búsqueda de dominios lo dividiremos en tres etapas:

- #### Envío de las secuencias
- #### Consulta del status del envío
- #### Descarga de resultados

### Cada etapa tiene una URL específica la cual definiremos a continuación
```python  
submit_url   = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
progress_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status"
results_url  = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result"
```  
### A cada etapa agregar headers específicos
```python  
submit_headers   = {"Accept":"text/plain"}
progress_headers = {"Accept":"text/plain"}
results_headers  = {"Accept":"text/tab-separated-values"}
```
## 7.1 Envío de las secuencias

### En esta etapa, construiremos un diccionario de python que adjuntaremos a `requests` para buscar los dominios funcionales
```python  
submit_data = {"email":"pamalemazac13@gmail.com",
               "title":"operon_0_05",
               "goterms":"false",
               "pathways":"false",
               "stype":"p",
               "sequence":query_str}
```
### Con los datos listos enviamos nuestra solicitud a InterProScan
  
submit_request = requests.post(submit_url,data=submit_data,headers=submit_headers)
  
### La API de InterProScan entrega un código de estado y un `job_id`.

> ### El código de salida del servidor web, nos indican si la solicitud fue exitosa:
>> * 1xx informational response – the request was received, continuing process
>>* 2xx successful – the request was successfully received, understood, and accepted
>> * 3xx redirection – further action needs to be taken in order to complete the request
>> * 4xx client error – the request contains bad syntax or cannot be fulfilled
>> * 5xx server error – the server failed to fulfil an apparently valid request

### Si el servidor nos entrega un código de estado 200, nuestras secuencias entraron en el servicio web

### Con el siguiente código podemos obtener tanto nuestro `job_id` como el código de estado
```python  
submit_status_code = submit_request.status_code
submit_job_id      = submit_request.text
print(submit_status_code)
print(submit_job_id)
  
200
iprscan5-R20230316-184520-0311-86547829-p2m
```  
## 7.2 Progreso del análisis de las secuencias
  
### Con nuestro `job_id` podemos consultar el progreso del análisis con las siguientes líneas de código
```python  
progress_request     = requests.get(progress_url+"/"+submit_job_id,headers=progress_headers)
progress_status_code = progress_request.status_code
progress_status      = progress_request.text
print(progress_status_code)
print(progress_status)
  
200
FINISHED
```  
## 7.3 Obtención de los resultados

### Si el estado del análisis es "FINISHED", podemos consultar el resultado de nuestro análisis

- ### `log` para acceder al reporte de texto del programa
- ### `tsv` para acceder al reporte tabular programa
```python  
results_log_request = requests.get(results_url+"/"+submit_job_id+"/log",headers=results_headers)
results_tsv_request = requests.get(results_url+"/"+submit_job_id+"/tsv",headers=results_headers)
  
print(results_log_request.text)
```  
![Captura de pantalla 2023-03-16 a la(s) 1 47 49 p m](https://user-images.githubusercontent.com/117956356/225736323-f783c407-c682-49cd-92ca-daf4db8695f3.png)

## 7.4 Examinación de los datos de dominios conservados

### Finalmente, usaremos la biblioteca StringIO para incorporar nuestros resultados de `requests` en un nuevo dataframe
```python  
results_tsv_str = StringIO(results_tsv_request.text)
results_column_names = ["sequence","md5","length","database","accession","description","start","end","evalue","post_processed","date","entry","name"]
results_df = pd.read_csv(results_tsv_str,sep="\t",names=results_column_names)
  
results_df
results_df["name"].unique()
```  
# **8. Conclusiones del ejercicio**
  
## En este ejercicio se lograron emplear correctamente bibliotecas que nos permitieron:

- ### Descargar secuencias genómicas de GenBank
- ### Predecir secuencias codificantes en el genoma de Klebsiella
- ### Obtener secuencias de referencia de UniProt asociadas a la resistencia a los antibióticos, el cual es un proceso molecular complejo que puede involucrar la transferencia de genes, la mutación y la selección natural
- ### Comparar las secuencias de referencia contra nuestras predicciones para determinar si en nuestro genoma hay genes codificantes asociados a la resistencia a los antibióticos que especificamos en la búsqueda de secuencias de UniProt
- ### Visualizar los datos genómicos a nivel grueso (pyCirclize) y a nivel fino (seaborn)
- ### Filtrar nuestros resultados para analizar unicamente un subset de secuencias, potencialmente asociadas a la resistencia a los antibióticos
- ### Búsqueda de dominios funcionales en el subset de proteínas de interés
  
  
  
  
  
  
  
  
  
  
  
