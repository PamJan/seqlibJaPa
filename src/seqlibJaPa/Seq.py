# definition of a class describing a sequence
class Seq:

    # dictonary witch complementary notation
    ALPH = {
       'A' : 'T',   'T' : 'A',   'G' : 'C',   'C' : 'G',
       'K' : 'M',   'M' : 'K',   'B' : 'V',   'D' : 'H',
       'H' : 'D',   'V' : 'B',   'N' : 'N',   '-' : '-'
    }
    
    # special method initialising a new object
    def __init__(self, seqid, title, seq):
        
        self.seqid = seqid
        self.title = title
        self.seq   = seq
        
    # special method returning the object length
    def __len__(self):
        return len(self.seq)
    
    # class method serialising sequences from a file
    @classmethod
    def from_file(cls, filename):

        # create an empty dictionary for seqid : Seq-object pairs
        # and a variable for the current seqid (None)
        seqs, seqid = {}, None

        # open the input file for reading in the text mode
        f = open(filename)

        # iterate over input file lines
        for line in f:
            
            # the first character being '>' means a new
            # sequence starts, extract seqid, title (if any),
            # create a new Seq object and place it
            # in the seqs dictionary
            if line[0] == '>':
                
                # [1:-1] -> skip the first '>' and the last '\n' (eol) character
                # split(' ', 1) -> split at the first occurrence of ' '
                header = line[1:-1].split(' ', 1)
                
                # seqid will always be the first element of the header list
                seqid  = header[0]
                
                # checking if seqid is unique
                if seqid in seqs:
                    raise Exception(f'Non-unique: {seqid}')
                    
                # if there is no space in the header, there is no title ('')
                # and the header list lengths is 1, instead of 2
                title  = header[1] if len(header) > 1 else ''
                
                # the seq property is initially set to an emtpy list
                seqs[seqid] = cls(seqid, title, [])
                
            else:
                # checking if character in IUPAC code
                if not set(line[:-1]).issubset(cls.ALPH.keys()):
                    raise Exception(f'character not in IUPAC code')
                # collect subsequent lines of a sequence from
                # as the seq list elements
                seqs[seqid].seq.append(line[:-1])

        # for each Seq object concatenate the lines
        # from the seq variable into one string variable
        for seq in seqs.values():
            seq.seq = ''.join(seq.seq)

        f.close()
            
        return seqs
        
    # The special method __repr__() returns a string representation of
    # an object, we will define this representation as 10 first letters
    # of the sequence stored within an object and three dots '...'
    # UNLESS the sequence length is shorter or equal to 10 letters.
    def __repr__(self):
        if(len(self.seq) <= 10):
            return self.seq
        else:
            return f'{self.seq[:10]}...'
            
    # special method, which returns a string whenever a DNASeq 
    # object is being converted to one,eg. str(obj) or print(obj).
    # it return FASTA formatted sequence string
    # sequence 60-character lines
    def __str__(self):
        # splitting string into lines
        lines = '\n'.join(
           [self.seq[i:i+60] for i in range(0, len(self.seq), 60)]
        )

        title = f' {self.title}' if self.title != '' else ''

        fasta = f'>{self.seqid} {title}\n{lines}'

        return fasta

    # method, witch return a new object of DNASeq type, which
    # will contain a reverse-complement sequence to the
    # one contained within the original object the method is called from.
    # '_revcmpl' suffix is added to the seqid from the original one
    def revcmpl(self):
        # revering the sequence
        seq = self.seq[::-1]
        
        # Using string method join(), the class dictionary ALPH and a
        # list comprehension expression, translating the reversed sequence and
        # converting into a string
        seq_revcmpl = ''.join([self.ALPH[el] for el in seq])
        
        # Creating seqid variable and assigning to it the object's seqid
        # together with the suffix '_revcmpl'.
        seqid = f'{self.seqid}_revcmpl'
        
        # Creating a new object of the DNASeq
        return DNASeq(seqid, self.title, seq_revcmpl)

    # special method __getitem__()
    # that allows to program what happens when an object
    # is indexed or sliced (like a list or a string),
    # eg. obj[3] or obj[4:5].
    # notation GeneBank!!!
    # - indexing starts at 1
    # - when a fragment is requested, both indices are inclusive:
    #   seq = 'ATGCTACG', seq[1:3] -> 'ATG'
    #          12345678
    # - the start index greater than stop index indicates
    #   a reverse complement (the complementary strand):
    #   seq[3:1] -> 'CAT'
    def __getitem__(self, key):
        
        if isinstance(key, slice):            
            start, end, step = key.start, key.stop, key.step
            
            if step is not None:
                # There is not an equivalent of step in GenBank notation
                
                raise KeyError('Step is not allowed in GenBank notation')
            
            if start is None:
                # Start must be provided.
                
                raise KeyError('start index is required')
            
            if end is None:
                # As well as the end.
                
                raise KeyError('end index is required')
                
            if not np.issubdtype(type(start), np.integer) or \
              not np.issubdtype(type(end), np.integer):
                # Start and end must be defined as integer values
                
                raise TypeError('Start and end must be integers')
            
            if start <= 0 or end <= 0:
                # Both must be greater than 0 as in GenBank notation
                # indexing starts at 1.
                
                raise KeyError('Minimal value for start and end is 1')
                

            # start and end are both integer type values and equal or greater to 1.

            if (start <= end):
                strand = 1
            else:
                strand = -1
                
            # If strand is equal to -1, swap the values of start and end
            if (strand == -1):
                a, b = (start, end)
                end = a
                start = b
            
            # new  seqid
            seqid = f'{self.seqid}_loc({start}_{end})'
            
            # Decreaseing start by 1. If start in GenBank notation is 1,
            # it is and equivalent of 0 in Python indexing (etc.),
            # so it needs to be decreased by 1 to be translated
            # from GenBank to Python.
            start -= 1
            
            # Create a new object of DNASeq witch new seqid 
            # and sliced
            sub_seq = DNASeq(seqid, self.title, self.seq[start:end])
            
            # If strand is -1, assign to the same reference ("variable")
            # a reverse complement of the new sequence object by using
            # the method revcmpl()
            
            if strand == -1:
                sub_seq = sub_seq.revcmpl()
                
            # returning the reference to the new sequence object. 
            return sub_seq
                
        else:
            
            # If we are here, it means that key is not a slice object.
            # Then we allow it be an integer greater than 0,
            # compliant with the GenBank notation.
            
            if not np.issubdtype(type(key), np.integer):
                raise TypeError('Index must be an integer')
                
            if key <= 0:
                raise KeyError('Minimal value of index is 1')
                
            # In case it is just one letter not a slice, we
            # will return simply a letter from the sequence
            # (a string value), NOT a new DNASeq object.
            # We need to remember to decrease the key value
            # by one to translate it from GenBank to Python.
            return self.seq[key-1]

        # The special method which will be invoked when two objects of the DNASeq 
        # type are being added, eg. seq_3 = seq_1 + seq_2. 
        def __add__(self, other):
            # seqid will be seqids of added objects,
            # separated by an underscore
            seqid = f'{self.seqid}_{other.seqid}'
            
            # title will be titles of added objects,
            # separated by an underscore
            title = f'{self.title}_{other.title}'

            # seq will simply be a concatenation of both sequences
            # stored within the added objects
            seq = f'{self.seq}{other.seq}'

            # returning the reference to the new object.
            return DNASeq(seqid, title, seq)

        # custom method copy() that will return an exact copy
        # of the existing object (as a type(self) type object).
        def copy(self):
            return DNASeq(self.seqid, self.title, self.seq)
            
# decide on what will be exposed as module content (here, only the Seq class)
__all__ = [ 'Seq' ]

