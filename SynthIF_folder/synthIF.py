#!/usr/bin/env python
# coding: utf-8

# In[18]:


import xml.etree.ElementTree as ET #parsování XML
import Bio #biopython
from Bio import SeqIO #parsování fasta souborů
from Bio.SeqRecord import SeqRecord
from Bio import Align #pairwise aligner
import re #regex
import subprocess #spouštění shell skriptů z prostředí pythonu
import os #přístup k funkcím operačního systému z pythonu
import random                                      
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from hmmlearn import hmm
import configparser


# In[19]:


#konfigurace
config = configparser.ConfigParser()
config.sections()
config.read('config.ini')
#geny
genes=[]
for gene, value in config['GENES'].items():
    if config['GENES'].getboolean(gene):
        genes.append(gene.upper())


# In[20]:


#GENEROVÁNÍ DAT - VYGENERUJE FASTA SOUBOR OBSAHUJÍCÍ ZVOLENÉ KIR GENY


output_file = "./outputTemp/geneSampleHT.fasta"

sequence = ""

for gene in genes:
    filename = f"sampleRefs/{gene}.fasta"

    for record in SeqIO.parse(filename, "fasta"):
        sequence += str(record.seq)

with open(output_file, "w") as f:
    f.write(">KIRpack\n")
    f.write(sequence + "\n")


# In[21]:


#GENEROVÁNÍ READŮ - PBSIM3

pbsim = subprocess.run('pbsim --strategy wgs --method errhmm --errhmm ~/miniconda3/data/ERRHMM-SEQUEL.model --genome ~/SynthIF_folder/outputTemp/geneSampleHT.fasta --depth 60 --length-mean 1429 --pass-num 10 --prefix ~/SynthIF_folder/outputTemp/reads', shell=True)


# In[22]:


#PŘEVOD NA FASTQ KVŮLI MINIMAP2 VSTUPU
samtools1 = subprocess.run('samtools fastq ~/SynthIF_folder/outputTemp/reads_0001.bam > ~/SynthIF_folder/outputTemp/reads.fastq', shell=True)


# In[23]:


#ZAROVNÁNÍ - MINIMAP2

output_file = "./outputTemp/mappingRef.fasta"

sequence = ""

for gene in genes:
    filename = f"mappingRefs/{gene}.fasta"

    for record in SeqIO.parse(filename, "fasta"):
        sequence += str(record.seq)

with open(output_file, "w") as f:
    f.write(">KIRpack\n")
    f.write(sequence + "\n")


minimap = subprocess.run('minimap2 -ax map-pb ~/SynthIF_folder/outputTemp/mappingRef.fasta ~/SynthIF_folder/outputTemp/reads.fastq > ~/SynthIF_folder/outputTemp/align.sam', shell=True)


# In[24]:


#SORTING
samtoolssort = subprocess.run('samtools sort ~/SynthIF_folder/outputTemp/align.sam -o ~/SynthIF_folder/outputTemp/align.sorted.bam', shell=True)


# In[25]:


#INDEXACE
samtoolsindex = subprocess.run('samtools index ~/SynthIF_folder/outputTemp/align.sorted.bam', shell=True)


# In[26]:


#KONSENZUÁLNÍ SEKVENCE - BCFTOOLS
bcf_pileup = subprocess.run('bcftools mpileup -Ou -f ~/SynthIF_folder/outputTemp/mappingRef.fasta ~/SynthIF_folder/outputTemp/align.sorted.bam > ~/SynthIF_folder/outputTemp/pileup.bcf', shell=True)
bcf_call = subprocess.run('bcftools call --ploidy 1 -mv -Oz -o ~/SynthIF_folder/outputTemp/calls.vcf.gz ~/SynthIF_folder/outputTemp/pileup.bcf', shell=True)
bcf_idx = subprocess.run('bcftools index ~/SynthIF_folder/outputTemp/calls.vcf.gz', shell=True)
bcf_cons = subprocess.run('bcftools consensus -f ~/SynthIF_folder/outputTemp/mappingRef.fasta ~/SynthIF_folder/outputTemp/calls.vcf.gz > ~/SynthIF_folder/outputTemp/consensus.fa', shell=True)


# In[27]:


#SIMULACE IDENTIFIKACE PŘÍTOMNÝCH GENŮ
lengths=[]
for gene in genes:
    filename = f"mappingRefs/{gene}.fasta"
    for record in SeqIO.parse(filename, "fasta"):  #zajištění délek referenčních sekvencí, vůči kterým bylo zarovnáno
        lengths.append(len(record.seq))

# načtení spojené sekvence
merged_record = next(SeqIO.parse("./outputTemp/consensus.fa", "fasta"))
sequence = merged_record.seq

# rozdělení
new_records = []
start = 0

for gene, length in zip(genes, lengths):
    end = start + length

    part_seq = sequence[start:end].upper()

    new_record = SeqRecord(
        part_seq,
        id=gene,
        description=""
    )

    new_records.append(new_record)
    start = end


SeqIO.write(new_records, "./outputTemp/split.fasta", "fasta")


# In[28]:


targets=[]
pattern = r"KIR\d+D[LS]\d+[AB]?"
targets_available=[]
for record in SeqIO.parse("./outputTemp/split.fasta", "fasta"): #načtení popsaného balíku KIR (fasta)
  targets.append(record.seq)
  match = re.search(pattern, record.description)
  if match:
    targets_available.append(match.group(0))


# In[29]:


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


# In[33]:


#TRÉNOVÁNÍ A SIMULACE HMM
def simul_hmm(target_idx, exony, n_components=8):
    gene=targets_available[target_idx]

    df=pd.read_csv("exonsBase/trSeq.csv")

    sekvence_train = (
        df[gene]
        .dropna()
        .astype(str)
        .tolist()
    )

    X = np.concatenate([ #převod na 0/1
        np.array([0 if c == "I" else 1 for c in seq])
        for seq in sekvence_train
    ]).reshape(-1, 1)

    lengths = [len(seq) for seq in sekvence_train]

    model = hmm.CategoricalHMM(  #diskrétní emise
        n_components=n_components,     # počet skrytých stavů
        n_iter=100,
        random_state=5
    )
    model.fit(X, lengths)
    pi = model.startprob_

    P = model.transmat_

    E = model.emissionprob_

    pocty = {}

    opakovani = 100000
    kroky = lengths[0]
    for pokus in range(opakovani):

        # výběr počátečního skrytého stavu podle pi
        r = random.random()
        soucet = 0

        for i in range(kroky):
            soucet += pi[i]
            if r < soucet:
                stav = i
                break

        vystup = ""

        for i in range(kroky):
            # emise
            r = random.random()  #náhodné číslo z uniformního rozdělení (0-1)
            if r < E[stav][0]:   #rozhodnutí mezi emitujícím stavem I/S - pokud je gen. číslo menší než ppst emise I, pak se emituje S
                vystup = vystup + "I"  #uloží I do sekvence stavů
            else:
                vystup = vystup + "S"  #jinak do sekvence uloží S

            # přechod do dalšího stavu
            r = random.random()    #náh. číslo
            soucet = 0    #proměnná pro hranici intervalu

            for novy_stav in range(8):
                soucet = soucet + P[stav][novy_stav]  #pravděpodobnost přechodu ze stavu stav do novy_stav - radek-sloupec matice P
                if r < soucet:  #testuji postupně všechny hranice pro vstup do daného stavu (transition FW, BW, self) - jiné stavy než propojené zajištěné nulou v matici přechodů
                    stav = novy_stav 
                    break #pokud se přejde do nového stavu, zastav a opakuj rozhodnutí o emisi

        if vystup in pocty:  #pokud je výstupní sekvence emitujících stavů v poli výstup, přičti počet
            pocty[vystup] = pocty[vystup] + 1
        else:  #pokud ne, pak zaveď nový a nastav počet na 1
            pocty[vystup] = 1

    #print("Počty jednotlivých posloupností:")  #výpisy
    serazene = sorted(pocty.items(), key=lambda x: x[1], reverse=True)

    for posloupnost, pocet in serazene:
        pravdepodobnost = pocet/opakovani
        print(posloupnost, ":", pocet, pravdepodobnost)

    rd=[posloupnost[0] for posloupnost in serazene]

    #print(rd)
    return rd, serazene, opakovani


# In[31]:


def export(gene, rd, serazene, exony, opakovani):
    output=[]

    for cd in rd:
        newIF=[]
        idx=0
        for idx, sym in enumerate(cd):
            if sym=="I":
                newIF.append(exony[idx])

        output.append(newIF)

    with open(f"outputs/isoforms{gene}.fasta", "w") as f:
        for i, outer in enumerate(output, start=1):
            m = "".join(str(seq) for seq in outer)
            pravdepodobnost = serazene[i-1][1]/opakovani

            f.write(f">{gene}_{i}_{pravdepodobnost}\n")
            f.write(m + "\n")


# In[32]:


results={}
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




# In[ ]:




