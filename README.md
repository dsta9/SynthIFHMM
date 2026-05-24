# SynthIFHMM
## Požadavky

### Operační systém
- Linux (testováno na Ubuntu 22.04 LTS nebo novější)

### Python prostředí
- Python 3.10 nebo novější
- Conda nebo Miniconda

### Použité Python knihovny
- Biopython
- NumPy
- Pandas
- Matplotlib
- hmmlearn

### Externí nástroje
- PBSim3
- SAMtools
- Minimap2
- BCFTools

Doporučená instalace bioinformatických nástrojů je pomocí Condy (`conda-forge` + `bioconda`).

---

# Instalace

## Klonování repozitáře

```bash
git clone https://github.com/dsta9/SynthIFHMM.git
cd SynthIFHMM
```

## Vytvoření prostředí

```bash
conda create -n synthifhmm python=3.10
conda activate synthifhmm
```

## Instalace Python knihoven

```bash
pip install biopython numpy pandas matplotlib hmmlearn
```

## Instalace bioinformatických nástrojů

```bash
conda install -c bioconda -c conda-forge pbsim3 samtools minimap2 bcftools
```

---

# Spuštění

Projekt se spouští pomocí:

```bash
python synthIF.py
```

V kořenovém adresáři musí být přítomen:
- `SynthIF.py`
- `config.ini`

---

# Konfigurace

Soubor `config.ini` slouží k nastavení:
- přítomnosti genů ve vstupních datech
- parametrů PBSIM3
- parametrů HMM

Hodnoty:
```ini
true  = gen bude přítomen
false = gen nebude přítomen
```

---

# Doporučení

Před spuštěním je doporučeno vyčistit pracovní adresáře:

```bash
./clearFolders.sh
```

> Skript smaže i obsah složky `outputs`.

---

# Parametry

| Parametr | Popis | Výchozí hodnota |
|---|---|---|
| `depth` | Průměrná hloubka pokrytí | `60` |
| `length` | Průměrná délka readu | `1429` |
| `hmmComponents` | Počet stavů HMM | `8` |
| `opakovani` | Počet simulací HMM | `100000` |
| `trIter` | Počet iterací trénování HMM | `100` |

---
