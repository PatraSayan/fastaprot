# fastaprot

A lightweight python library for extracting useful information from uniprot fasta files.

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/PatraSayan/fastaprot.git
cd fastaprot
pip install -e .
```

It installs `numpy` and `matplotlib` automatically as dependencies.

## Usage

### Reading a UniProt FASTA file

```python
from fastaprot import read_fasta

data = read_fasta("tests/INSR.fasta")
print(data["header"])
print(data["sequence"])
```
### Extracting metadata

```python
from fastaprot import meta_info

info = meta_info(data["header"])
print(info)
```

### Sequence Information Extraction

```python
from fastaprot import Sequence_Info

seq_info = Sequence_Info(data["sequence"])

# Get length of the Protein
print(seq_info.length())

# Get number of different amino acids in the sequence
print(seq_info.amino_acid_count())

# Plot percentage of different amino acids
print(seq_info.amino_acid_plot())

# Get a theoretical molecular weight 
print(seq_info.molecular_weight())

# Get number of negative charged amino acid (Asp + Glu)
print(seq_info.neg_amino_acid())

# Get number of positive charged amino acid (His + Lys + Arg)
print(seq_info.pos_amino_acid())

# Get tryptophan count
print(seq_info.tryptophan())
```

### Enzymatic Digestion

```python
from fastaprot import Digest

digest = Digest(data["sequence"])

# trypsin digestion (cleaves after K/R, skips K/R-P bonds)
print(digest.trypsin())  

# trypsin/p digestion (cleaves after K/R)
print(digest.trypsin_p())

# chymotrypsin digestion (cleaves after F/W/Y)
print(digest.chymotrypsin())

# CNBr digestion (cleaves after M)
print(digest.cnbr())
```

## Project Structure

```
── src/fastaprot/
│   ├── __init__.py
│   └── main.py
├── tests/
│   ├── INSR.fasta
│   └── test.py
├── LICENSE
├── pyproject.toml
└── README.md
```

## Requirements
- Python >= 3.9
- numpy
- matplotlib

## License
- MIT License

## Author
Sayan Patra ([patrasayan2024@gmail.com](mailto:patrasayan2024@gmail.com))