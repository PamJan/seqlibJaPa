# The library overview

seqlibJaPa is a Python library for reading and analyzing fasta files witch GeneBank notatnion. 
It can be used for:
- reading multiple sequnces from fasta files
- geting reverse-complement sequence
- slicing and indexing sequence in GeneBank notation - list[start:stop:step]
- concatenating and coping sequences

# Library installation using _pip_
Installation of the `seqlibJaPa` library with pip is quite straightforward:
```Bash
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple seqlibJaPa
```

# A quick example of the library usage
Learn how to use the `seqlibFiLa` library with examples provided below.

```Python
# Import the Seq class from the seqlibJaPa library
from seqlibJaPa import Seq

# Example 1 - creating DNASeq object from scratch
seq = DNASeq('new_seq', 'an exemplary sequence',
             'ATCGTAGGATCGGATTAGAGCGATTAGCTAG')

print(f'{seq.seqid} ({seq.title}): {seq.seq[:5]}...')

# Example 2 - reading fasta file

# Reload the sequences to have a collection of objects
# that are instances of the up-to-date DNASeq class.

seqs = DNASeq.from_file('input/Staphylococcus_MLST_genes.fasta')

seqs


# Example 3 - lenght of selected sequence

# Reload the sequences to have a collection of objects
# that are instances of the up-to-date DNASeq class:

seqs = DNASeq.from_file('input/Staphylococcus_MLST_genes.fasta')
# Select one of the sequences by its sequence id (seqid):
seq = seqs['yqiL']

# Look up the length of the contained sequence:
len(seq)


# Example 4 - slicing

# Reload the sequences to have a collection of objects
# that are instances of the up-to-date DNASeq class:

seqs = DNASeq.from_file('input/Staphylococcus_MLST_genes.fasta')

# Select 10 nucleotide fragments of two different
# sequences and add them together:
seq_1 = seqs['arcC'][1:10]

seq_2 = seqs['glpF'][1:10]

seq_3 = seq_1 + seq_2

# Print the initial sequences and the result of their addition
# as if they were string values:
print(seq_1, seq_2, seq_3, sep='\n\n')

```

# The `DNASeq` class in details
DNASeq is a class that storage: 
- `seqid` - ID of the sequence eg. "yqiL"
- `title` - sequence title eg. "Acetyle coenzyme A acetyltransferase"
- `seq` - sequence eg "TTGCTGTGCACGTACTGCTTTTTGTT..."

It is posible to initialise it by monually adding sequence or by reading sequence/s from fasta files 
 <br /> -> look **[Initialisation](#Initialisation)**

DNAseq class have one variable `ALPH` witch is a dictonary witch complementary notation. <br /><br />
It have implemented 9 methods
 <br /> -> look **[Methods](#Methods)**

```Python
class DNASeq:

    ALPH
    
    __init__()
    __repr__()
    __str__()
    __len__()
    __add__()
    __getitem__()
    
    from_file()
    revcmpl()
    copy()
```


## Initialisation
Initialising one sequence
```Python
seq = DNASeq('new_seq', 'an exemplary sequence',
             'ATCGTAGGATCGGATTAGAGCGATTAGCTAG')
```            

Reding sequences from fasta file
```Python
seq = DNASeq('new_seq', 'an exemplary sequence',
             'ATCGTAGGATCGGATTAGAGCGATTAGCTAG')
```
## Methods
- `__init__(self, seqid, title, seq)` &ndash; Initialisation. Returns `DNASeq` object. <br />  -> look **[Initialisation](#Initialisation)**
- `__repr__(self)` &ndash; The special method which returns a string representation of an object. Returns 0 first letters of the sequence stored within an object and three dots '...'UNLESS the sequence length is shorter or equal to 10 letters.
- `__str__(self)` &ndash; The special method which returns a string whenever a DNASeq object is being converted to one, eg. str(obj) or print(obj). Returns the sequence contained within an object as FASTA formatted sequence string.
- `__len__(self)` &ndash; The special method which returns the length of the sequence contained within a given DNASeq object when the reference to that object is passed to the built-in len() function. Returns lenght of sequence.
- `__getitem__(self, key)` &ndash; The special method which will be invoked when object is indexed or sliced (like a list or a string),eg. obj[3] or obj[4:5] in GeneBank notation. <br />  -> look **[Indexing and slicing](#Indexing-and-slicing)**
- `__add__()` &ndash; The special method which will be invoked when two objects of the DNASeq type are being added, eg. seq_3 = seq_1 + seq_2. Returns reference to the new object <br />  -> look **[Adding two sequences](#Adding-two-sequences)**
- `from_file(self, filename)` &ndash; Deserialisation of sequences from a FASTA format file to a collection of DNASeq objects. Returns `DNASeq` objects.
- `revcmpl(self)` &ndash; Method creating reverse-complement sequence. Returns a new object of DNASeq type, which will contain a reverse-complement sequence to the one contained within the original object the method is called from. '_revcmpl' suffix is added to the seqid from the original one.
- `copy(self)` &ndash; Metchod return an exact copy of the existing object (as a type(self) type object).

## Indexing and slicing
`__getitem__(self, key)` &ndash; The special method which will be invoked when object is indexed or sliced (like a list or a string),eg. obj[3] or obj[4:5] in GeneBank notation:
- indexing starts at 0
- when a fragment is requested, both indices are inclusive:  seq = 'ATGCTACG', seq[1:3] -> 'ATG'
- the start index greater than stop index indicates a reverse complement (the complementary strand): seq[3:1] -> 'CAT'
Returns sliced sequence.

## Adding two sequences
`__add__(self)` &ndash; The special method which will be invoked when two objects of the DNASeq type are being added, eg. seq_3 = seq_1 + seq_2. Returns reference to the new object where:
- `seqid` is seqids of added objects, separated by an underscore;
- `title` is titles of added objects separated by an underscore;
- `seq` is simply be a concatenation of both sequences stored within the added objects.

