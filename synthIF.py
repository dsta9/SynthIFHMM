#!/usr/bin/env python
# coding: utf-8

# In[2]:


import xml.etree.ElementTree as ET #parse XML
import Bio #biopython
from Bio import SeqIO #parsing fasta souborů
from Bio.SeqRecord import SeqRecord #práce s FASTA hlavičkami
from Bio import Align #pairwise aligner
import re #regex
import subprocess #spouštění shell skriptů z prostředí pythonu
import os #přístup k funkcím operačního systému z pythonu
import random  #RNG
import matplotlib.pyplot as plt # grafy
import numpy as np #matematika
import pandas as pd #tabulky
from hmmlearn import hmm #HMM
import configparser #konfigurační soubory


# In[3]:


# načtení konfigurace
config = configparser.ConfigParser()
config.sections()
config.read('config.ini')
#načtení zvolených genů
genes=[]
for gene, value in config['GENES'].items():
    if config['GENES'].getboolean(gene):
        genes.append(gene.upper())


# In[4]:


#GENEROVÁNÍ DAT - VYGENERUJE FASTA SOUBOR OBSAHUJÍCÍ ZVOLENÉ KIR GENY


output_file = "./outputTemp/geneSampleHT.fasta" #výstupní soubor

sequence = ""  #výstupy k exportu

for gene in genes:    #pro každý zvolený gen zařadí jeho sekvenci do balíku
    filename = f"sampleRefs/{gene}.fasta"

    for record in SeqIO.parse(filename, "fasta"):
        sequence += str(record.seq)

with open(output_file, "w") as f:
    f.write(">KIRpack\n")   #sekvnce se dočasně jmenuje KIRpack
    f.write(sequence + "\n")


# In[60]:


#GENEROVÁNÍ READŮ - PBSIM3

#parametry generátoru
depth = config['PBSIM']['depth'] #průměrná hloubka pokrytí
length = config['PBSIM']['length'] #průměrná délka readu

#spuštění generátoru modelbased sim, HMM model chyby, params
pbsim = subprocess.run('pbsim --strategy wgs --method errhmm --errhmm ~/miniconda3/data/ERRHMM-SEQUEL.model --genome ~/SynthIFHMM/outputTemp/geneSampleHT.fasta --depth '+str(depth)+' --length-mean '+str(length)+' --pass-num 10 --prefix ~/SynthIFHMM/outputTemp/reads', shell=True)


# In[61]:


#PŘEVOD NA FASTQ KVŮLI MINIMAP2 VSTUPU
samtools1 = subprocess.run('samtools fastq ~/SynthIFHMM/outputTemp/reads_0001.bam > ~/SynthIFHMM/outputTemp/reads.fastq', shell=True)


# In[62]:


#ZAROVNÁNÍ - MINIMAP2

output_file = "./outputTemp/mappingRef.fasta"

sequence = ""

#na základě apriorní znalosti zvolených genů se z referenčních sekvencí sestaví předloha pro zarovnání

for gene in genes:
    filename = f"mappingRefs/{gene}.fasta"

    for record in SeqIO.parse(filename, "fasta"):
        sequence += str(record.seq)

with open(output_file, "w") as f:
    f.write(">KIRpack\n")
    f.write(sequence + "\n")

#spuštění MINIMAP s nastavenými parametry pro PacBio

minimap = subprocess.run('minimap2 -ax map-pb ~/SynthIFHMM/outputTemp/mappingRef.fasta ~/SynthIFHMM/outputTemp/reads.fastq > ~/SynthIFHMM/outputTemp/align.sam', shell=True)


# In[63]:


#SORTING
samtoolssort = subprocess.run('samtools sort ~/SynthIFHMM/outputTemp/align.sam -o ~/SynthIFHMM/outputTemp/align.sorted.bam', shell=True)


# In[64]:


#INDEXACE
samtoolsindex = subprocess.run('samtools index ~/SynthIFHMM/outputTemp/align.sorted.bam', shell=True)


# In[65]:


#KONSENZUÁLNÍ SEKVENCE - BCFTOOLS
bcf_pileup = subprocess.run('bcftools mpileup -Ou -f ~/SynthIFHMM/outputTemp/mappingRef.fasta ~/SynthIFHMM/outputTemp/align.sorted.bam > ~/SynthIFHMM/outputTemp/pileup.bcf', shell=True)
bcf_call = subprocess.run('bcftools call --ploidy 1 -mv -Oz -o ~/SynthIFHMM/outputTemp/calls.vcf.gz ~/SynthIFHMM/outputTemp/pileup.bcf', shell=True)
bcf_idx = subprocess.run('bcftools index ~/SynthIFHMM/outputTemp/calls.vcf.gz', shell=True)
bcf_cons = subprocess.run('bcftools consensus -f ~/SynthIFHMM/outputTemp/mappingRef.fasta ~/SynthIFHMM/outputTemp/calls.vcf.gz > ~/SynthIFHMM/outputTemp/consensus.fa', shell=True)


# In[5]:


#SIMULACE IDENTIFIKACE PŘÍTOMNÝCH GENŮ
lengths=[]
for gene in genes:
    filename = f"mappingRefs/{gene}.fasta"
    for record in SeqIO.parse(filename, "fasta"):  #zajištění délek referenčních sekvencí, vůči kterým bylo zarovnáno (protože apriorní info je k dispozici)
        lengths.append(len(record.seq))

# načtení spojené sekvence
merged_record = next(SeqIO.parse("./outputTemp/consensus.fa", "fasta"))
sequence = merged_record.seq

# rozdělení
new_records = []
start = 0

for gene, length in zip(genes, lengths):
    end = start + length

    part_seq = sequence[start:end].upper() #zjištěný úsek se musí převést na velká písmena

    new_record = SeqRecord(   #knihovna SeqRecord vytvoří genomický fasta soubor pomocí konstruktoru
        part_seq,
        id=gene,
        description=""
    )

    new_records.append(new_record)   #přidání záznamu do souboru již "rozpoznaných" KIR
    start = end


SeqIO.write(new_records, "./outputTemp/split.fasta", "fasta") #export do souboru SPLIT


# In[6]:


targets=[]
pattern = r"KIR\d+D[LS]\d+[AB]?" #vyhledá se konkrétní patern začínající KIR, následuje číslo, písmeno D, pak L, nebo S, dále může následovat A, nebo B
targets_available=[] #pole dostupných genů v sekvenci
for record in SeqIO.parse("./outputTemp/split.fasta", "fasta"): #načtení popsaného balíku KIR (fasta)
  targets.append(record.seq)
  match = re.search(pattern, record.description)
  if match:
    targets_available.append(match.group(0))  #ukládání do výstupního pole validních cílů


# In[7]:


#INPUTS: index genu uloženého v poli validních cílů
#OUTS: použitá reference, pole získaných exonů, best (výstup z Pairwise aligner), transkript_id (název genu, kterému patří exony) 

def extract_exons(target_idx):
    target_test=targets[target_idx] #výběr jednoho KIR z balíku (sekvence)
    transcript_id=targets_available[target_idx]

    db=ET.parse("./exonsBase/kir_base.xml") #parsování XML - výběr z databáze referencí exonů
    root=db.getroot()

    tc=root.find(f".//transcript[@id='{transcript_id}']")

    exons = []
    exon_elements=tc.findall("exon")
    for exon in exon_elements:
      seq=(exon.text or "").strip()
      exons.append(seq)
    reference=exons
    aligner=Align.PairwiseAligner()       #zarovnání sekvencí
    aligner.mode="local"
    best=[]
    aligner.match_score = 2
    aligner.mismatch_score = -2
    aligner.open_gap_score = -20    #penalizace pro odstranění nepřítomných exonů
    aligner.extend_gap_score = -5


    for i in range (0,len(reference)):   #postupné zarovnání známých exonů na vstupní sekvenci
      alignments = aligner.align(reference[i], target_test)
      if(alignments[0].score > 20):  #ošetření nepřítomných exonů -> skóre musí být větší než 30
        best.append(alignments[0])

    exony = []

    exony.append(
        target_test[best[0].coordinates[1][0]:best[0].coordinates[1][-1]]
    )

    for i in range(1, len(best)):
        start = best[i].coordinates[1][0]
        end = best[i].coordinates[1][1]

        exony.append(
            target_test[start:end]
        )

    """
    print(1)
    print(target_test[best[0].coordinates[1][0]:best[0].coordinates[1][-1]])  #vstupní sekvence se vezme od začátku až do referenčního konce exonu
    for i in range(1,len(best)):
      print(i+1)
      print(target_test[best[i].coordinates[1][0]:best[i].coordinates[1][-1]])  #rozložená vstupní sekvence na jednotlivé exony
      print(best[i].score)
    """
    return reference, exony, best, transcript_id


# In[8]:


#TRÉNOVÁNÍ A SIMULACE HMM
#INPUTS: index genu v poli validních cílů, pole exonů předané z metody extract_exons, počet komponent HMM (z konfigu, nebo není-li předáno, nastavuje se 8)
#OUTPUTS: rd (pole unikátních generovaných posloupností), serazene (pole dvojic posloupnost tokenů/absolutní četnost seřazené sestupně), opakovani (počet simulací)



def simul_hmm(target_idx, exony, n_components=8):
    gene = targets_available[target_idx] #nastavení aktuálně zpracovávaného záznamu

    df = pd.read_csv("exonsBase/trSeq.csv") #načtení trénovacích množin

    sekvence_train = (   
        df[gene]
        .dropna()
        .astype(str)
        .tolist()
    )

    # převod I/S na 0/1
    X = np.concatenate([
        np.array([0 if c == "I" else 1 for c in seq])
        for seq in sekvence_train
    ]).reshape(-1, 1)

    lengths = [len(seq) for seq in sekvence_train]

    model = hmm.CategoricalHMM(  #vytvoření modelu pomocí hmmLearn
        n_components=n_components,
        n_iter=config.getint('PARAMETERS', 'trIter'),   #počet iterací se načte z konfigu
        random_state=5    #seed
    )

    model.fit(X, lengths)   #trénování modelu (reestimace parametrů A, B, pi)

    pocty = {} #četnosti výstupů

    opakovani = config.getint('PARAMETERS', 'opakovani')  #počet simulací
    kroky = lengths[0] #počet časových okamžiků HMM

    for pokus in range(opakovani):
        # simulace pomocí hmmlearn
        X_gen, Z_gen = model.sample(kroky)

        # převod 0/1 na tokeny I/S
        vystup = "".join("I" if x[0] == 0 else "S" for x in X_gen)

        pocty[vystup] = pocty.get(vystup, 0) + 1

    serazene = sorted(pocty.items(), key=lambda x: x[1], reverse=True) #řazení podle četnosti sestupně

    for posloupnost, pocet in serazene:
        pravdepodobnost = pocet / opakovani
        print(posloupnost, ":", pocet, pravdepodobnost)

    rd = [posloupnost[0] for posloupnost in serazene]

    return rd, serazene, opakovani


# In[9]:


#INPUTS: gene (zpracovávaný gen), rd, serazene, exony, opakovani (viz výše)
#OUTPUTS: finální výstup v souboru fasta 

def export(gene, rd, serazene, exony, opakovani):
    output=[] #proměnná pro uložení výstupních informací

    for cd in rd:
        newIF=[]
        idx=0
        for idx, sym in enumerate(cd):
            if sym=="I":
                newIF.append(exony[idx])   #substituce tokenů I za příslušný exon daného genu v závislosti na jeho pořadí, pokud je token S, pak se přeskočí

        output.append(newIF)  #řetězení exonů = vznik izoformy

    with open(f"outputs/isoforms{gene}.fasta", "w") as f:
        for i, outer in enumerate(output, start=1):
            m = "".join(str(seq) for seq in outer)
            pravdepodobnost = serazene[i-1][1]/opakovani

            f.write(f">{gene}_{i}_{pravdepodobnost}\n") #zápis do fasta souboru
            f.write(m + "\n")


# In[10]:


results={}
#pro každý gen dle konfigurace se volá metoda extract_exons, předání výstupů do metody simul_hmm, uložení do prom. results, která se předá jako parametr metodě export
for target_idx in range(len(targets)):
    reference, exony, best, gene = extract_exons(target_idx)
    rd, serazene, opakovani = simul_hmm(target_idx, exony, config.getint('PARAMETERS', 'hmmComponents'))

    results[gene] = {
        "rd": rd,
        "serazene": serazene,
        "opakovani": opakovani,
        "exony": exony
    }


    export(gene, rd, serazene, exony, opakovani)



