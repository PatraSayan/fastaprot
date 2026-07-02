# essential library
import numpy as np
import matplotlib.pyplot as plt
import re
import os

# creating a output folder
os.makedirs("output", exist_ok=True)

# Reading the fasta and extracting the header and the sequence
def read_fasta(dir):
    sequence = ""
    with open(dir) as f:
        for line in f:
            line = line.strip()                     # strip to remove all white space

            if line.startswith(">"):
                header = line[1:]                   # taking the  header information
            else: 
                sequence += line                    # Actual sequence
    return {"header":header, "sequence": sequence}  # returnsing a dicsionary

# Getting additional information like organism, gene, version, exsistence
def meta_info(header):
    db, acession, rest = header.split("|",2)
    entry, _, rest = rest.partition(" ")
    additional = {}
    additional_pattern = re.compile(r"(OS|OX|GN|PE|SV)=(.*?)(?=\s+(?:OS|OX|GN|PE|SV)=|$)")

    for match in additional_pattern.finditer(rest):
        additional[match.group(1)] = match.group(2).strip()
    
    description = re.search(r"\s(?:OS|OX|GN|PE|SV)=", rest)
    first_index = description.start()
    description = rest[:first_index].strip()

    return {
        "Database": db,
        "Accession": acession,
        "Uniprot Entry": entry,
        "Organism": additional.get("OS"),
        "Taxonomy ID": additional.get("OX"),
        "Gene Name": additional.get("GN"),
        "Protein Existence": additional.get("PE"),
        "Sequence Version": additional.get("SV")
    }

# Amino acid single letter code
amino_acids = {
    "Alanine": "A",
    "Arginine": "R",
    "Asparagine": "N",
    "Aspartic acid": "D",
    "Cysteine": "C",
    "Glutamic acid": "E",
    "Glutamine": "Q",
    "Glycine": "G",
    "Histidine": "H",
    "Isoleucine": "I",
    "Leucine": "L",
    "Lysine": "K",
    "Methionine": "M",
    "Phenylalanine": "F",
    "Proline": "P",
    "Serine": "S",
    "Threonine": "T",
    "Tryptophan": "W",
    "Tyrosine": "Y",
    "Valine": "V"
}

# estimated molecular weight 
aa_mw = {
    "A": 89.09318,
    "R": 174.2017,
    "N": 132.1179,
    "D": 133.1027,
    "C": 121.1590,
    "Q": 146.1445,
    "E": 147.1293,
    "G": 75.0666,
    "H": 155.1546,
    "I": 131.1736,
    "L": 131.1736,
    "K": 146.1882,
    "M": 149.2124,
    "F": 165.1900,
    "P": 115.1310,
    "S": 105.0930,
    "T": 119.1197,
    "W": 204.2262,
    "Y": 181.1894,
    "V": 117.1469
}

# digestion information
class Digest:

    def __init__(self, sequence):
        self.sequence = sequence
    
    # Trypsin/p total cleavage
    def trypsin_p(self):
        peptide = []
        start = 0

        for i, aa in enumerate(self.sequence):
            if aa in ("K","R"):
                peptide.append(self.sequence[start:i+1])
                start = i + 1
        if start < len(self.sequence):
            peptide.append(self.sequence[start:])
        
        return "Complete Trypsin/P digested peptides:", peptide
    
    def cnbr(self):
        peptide = []
        start = 0

        for i, aa in enumerate(self.sequence):
            if aa == "M":
                peptide.append(self.sequence[start:i+1])
                start = i+1
        if start < len(self.sequence):
            peptide.append(self.sequence[start:])
        
        return "Complete CNBr digested peptides:",peptide
    
    def chymotrypsin(self):
        peptide = []
        start = 0

        for i, aa in enumerate(self.sequence):
            if aa in ("F", "W", "Y"):
                peptide.append(self.sequence[start:i+1])
                start = i+1
        if start < len(self.sequence):
            peptide.append(self.sequence[start:])

        return "Complete chymotrypsin digested peptides:",peptide
    
    def trypsin(self):
        peptide = []
        start = 0

        for i, aa in enumerate(self.sequence):
            if aa in ("K","R") and (i==len(self.sequence) - 1 or self.sequence[i+1] != "P"):
                peptide.append(self.sequence[start:i+1])
                start = i + 1
        if start < len(self.sequence):
            peptide.append(self.sequence[start:])
        
        return "Complete trypsin digested peptides:",peptide



# Other information about the sequence
class Sequence_Info:

    def __init__(self, sequence):
        self.sequence = sequence

    def length(self):
        return len(self.sequence)
    
    def amino_acid_count(self):
        aa_list = {}
        for key, value in amino_acids.items():
            amino_acid = key
            aa_count = self.sequence.count(value)
            aa_percent = (self.sequence.count(value) / len(self.sequence)) *100
            aa_list[amino_acid] = (aa_count, aa_percent)
        lines = [f"{'Amino Acid':<15}{'Count':<10}{'Percentage':<10}", "-" * 35]
        for key in aa_list:
            lines.append(f"{key:<15}{aa_list[key][0]:<10}{aa_list[key][1]:<10.2f}")
        return "\n".join(lines)
    
    def amino_acid_plot(self):
        aa = {}
        for key, value in amino_acids.items():
            aa[key] = self.sequence.count(value)/len(self.sequence) * 100
        aa_asc = dict(sorted(aa.items(), key = lambda item:item[1]))
        plt.figure(figsize = (12,6))
        colours = plt.cm.Spectral(np.linspace(0,1,len(aa_asc)))
        plt.bar(aa_asc.keys(), aa_asc.values(), color = colours)
        plt.xticks(rotation = 45, ha = "right", fontweight="bold");
        plt.yticks(fontweight = "bold", fontsize = 14);
        plt.xlabel("Amino Acids", fontweight = "bold", fontsize = 14)
        plt.ylabel("Percentage in Fasta sequence", fontweight = "bold", fontsize = 14)
        plt.savefig("output/amino_acid_plot.png", dpi = 300, bbox_inches = "tight")
        plt.show()
    
    def molecular_weight(self):
        MW = 0
        for i in self.sequence:
            MW += aa_mw[i]
        return MW - (18.01528 * (len(self.sequence)-1))
    
    def neg_amino_acid(self):
        Asp = self.sequence.count("D")
        Glu = self.sequence.count("E")
        return Asp + Glu
    
    def pos_amino_acid(self):
        His = self.sequence.count("H")
        Lys = self.sequence.count("K")
        Arg = self.sequence.count("R")
        return Lys + Arg + His
    
    def tryptophan(self):
        return self.sequence.count("W")
    