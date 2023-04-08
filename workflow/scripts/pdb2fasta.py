
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

template_txt = str(snakemake.input)
fasta_file = str(snakemake.output)

aminoacids = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "LYR": "K"
}

def to_fasta(fname):
    has_seqres = False
    seq = []
    with open(fname) as fh:
        for line in fh:
            if line.startswith('SEQRES'):
                has_seqres = True
                tag, row_num, chain, length, *residues = line.split()
                if chain == 'A':
                    seq += [ aminoacids[x] for x in residues ]
            elif not has_seqres and line.startswith('ATOM'):
                tag, atom_num, atom, residue, chain, residue_num, *rest = line.split()
                if chain == 'A' and atom == 'N':
                    seq.append(aminoacids[residue])
    return ''.join(seq)

with open(fasta_file, 'w') as out:
    with open(template_txt) as fh:
        for line in fh:
            name, tag, pdb_file = line.split()
            seq = to_fasta(pdb_file)
            record = SeqRecord(Seq(seq), name[1:], description = '')
            SeqIO.write(record, out, 'fasta')
