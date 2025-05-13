'''
seqlibJaPa

DNASeq is a Python library for reading and analyzing fasta files witch GeneBank notatnion. 
It can be used for:
- reading multiple sequnces from fasta files
- geting reverse-complement sequence
- slicing and indexing sequence in GeneBank notation - list[start:stop:step]
- concatenating and coping sequences

'''

from seqlibFiLa.Seq import Seq

__version__ = '0.1.0'
__license__ = 'GPL-3.0-or-later'
__status__  = 'beta'
__author__  = 'Jan Pamuła'
__email__   = 'jan.pamula@student.uj.edu.pl'
__credits__ = [ 'Jagiellonian University, Krakow, Poland', 'PP4B Team' ]

__all__ = [
    '__version__',
    '__license__',
    '__status__',
    '__author__',
    '__email__',
    '__creadits__',
    'Seq'
]

