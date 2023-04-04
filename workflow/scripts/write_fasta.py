
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

with open(str(snakemake.output), 'w') as file:
    for name, seq in snakemake.params['seq'].items():
        record = SeqRecord(Seq(seq), id = name)
        SeqIO.write(record, file, 'fasta')
