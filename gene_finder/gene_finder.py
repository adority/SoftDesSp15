# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: Aubrey Dority

"""

# you may find it useful to import these variables (although you are not required to use them)
from amino_acids import aa, codons, aa_table
import random
from load import load_seq

def shuffle_string(s):
    """ Shuffles the characters in the input string
        NOTE: this is a helper function, you do not have to modify this in any way """
    return ''.join(random.sample(s,len(s)))

### YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
        
        I included more unit testing because it's important that the code find the right complement for each nucleotide and since there are only 4, it's easy enough to test that.

    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this
    if nucleotide == 'T':
        return 'A'
    elif nucleotide == 'A':
        return 'T'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'C':
        return 'G'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
    
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

        These unit tests are probably sufficient because they use strings of various lengths and letter combinations.

    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    dna2=''
    for nucleotide in dna[::-1]:
        dna2=dna2+get_complement(nucleotide)
    return dna2


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string

        I included another unit test because I wanted to be sure the whole string would be returned if there weren't a stop codon.

    >>> rest_of_ORF("ATGTACTTG")
    'ATGTACTTG'
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    i = 0
    stop_codons=['TAG','TGA','TAA']
    while i < len(dna):
        if (dna[i:(i+3)] in stop_codons):
            return dna[0:i]
        else:
            i=i+3
    return dna

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs

        This seems like a sufficient unit test since it includes a nested ORF in the original string and also stop codons that aren't in the frame of the sequence.

    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    i=0
    stop_codons=['TAG','TGA','TAA']
    ORFs=[]
    while i<len(dna):
        if dna[i:(i+3)] == 'ATG':
            start=i
            while not(dna[i:(i+3)] in stop_codons) and i<len(dna):
                i=i+3
            ORFs.append(dna[start:i])
        else:
            i=i+3
    return ORFs

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs

	This is a good unit test because it includes sequences in each frame.

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    ORFs=[]
    for i in range(3):
        ORFs+=find_all_ORFs_oneframe(dna[i:])
    return ORFs

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs

        This seems like a sufficient test as long as it produces results from both strands given one input string. 

    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    results=find_all_ORFs(dna)
    strand2=get_reverse_complement(dna)
    results2=find_all_ORFs(strand2)
    return results+results2

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string

        This is a sufficient test because since the longest ORF is on the other strand, we can assume it checked both.
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("TTTGATGCTACATTCGCAT")
    'ATGCTACATTCGCAT'
    """
    allORFs=find_all_ORFs_both_strands(dna)
    length1=1
    longestORF=''
    for ORF in allORFs:
        if len(ORF)>length1:
            length1=len(ORF)
            longestORF=ORF

    return longestORF


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 

        Seems like you can't really unit test this because the shuffles are random, but the lengths returned should be less than or equal to the length of the original string,
        approaching that length as num_trials increases.
        """

    i=0 
    length=1   
    while i<num_trials:
        dna=shuffle_string(dna)
        ORF=longest_ORF(dna)
        if len(ORF)>length:
            length=len(ORF)
        i+=1
    return length


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        I don't think any more unit tests are necessary because these include a string of a length that isn't a multiple of 3 as well as one that is,
        meaning there aren't problems with how the code deals with incomplete codons.

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    i=0
    a_as=''
    while i < (len(dna) - 2):
        codon = dna[i:(i + 3)]
        amino_acid = aa_table[codon]
        a_as += amino_acid
        i += 3
    return a_as

def gene_finder(dna):
    """ Returns the amino acid sequences coded by all genes that have an ORF
        larger than the specified threshold.
        
        dna: a DNA sequence
        threshold: the minimum length of the ORF for it to be considered a valid
                   gene.
        returns: a list of all amino acid sequences whose ORFs meet the minimum
                 length specified.
        Unit testing this will be difficult since longest_ORF_noncoding is going to return a threshold equal to the length of the initial dna strand for short strands;
        therefore none of the ORFs in the DNA will exceed the threshold. For longer strands of DNA, it's easier just to try the code and see if it returns reasonable strings
        of amino acids.
    """
    aminos=[]
    threshold=longest_ORF_noncoding(dna, 1500)
    ORFs=find_all_ORFs_both_strands(dna)
    for ORF in ORFs:
        if len(ORF)>threshold:
            aminos.append(coding_strand_to_AA(ORF))
    return aminos

sequence=load_seq("./data/X73525.fa")
gene_finder(sequence)

# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 22:02:04 2014

@author: pruvolo
"""

from os import path

def load_seq(fasta_file):
    """ Reads a FASTA file and returns the DNA sequence as a string.

    fasta_file: the path to the FASTA file containing the DNA sequence
    returns: the DNA sequence as a string
    """
    retval = ""
    f = open(fasta_file)
    lines = f.readlines()
    for l in lines[2:]:
        retval += l[0:-1]
    f.close()
    return retval
    
def load_salmonella_genome():
    f = open(path.join('.','data','salmonella_all_proteins'))
    lines = f.readlines()
    retval = []
    gene = []
    is_amino_acid_seq = False
    
    for line in lines:
        if line[5:].find("CDS") == 0:
            coords = line[21:-1]
            if len(gene) == 3:
                retval.append(gene)
            gene = [coords]
        elif line[21:].find("/protein_id") == 0:
            gene.append(line[34:-2])
        elif line[21:].find("/translation") == 0:
            if line[-2] != '"':
                amino_acid_seq = line[35:-1]
                is_amino_acid_seq = True
            else:
                amino_acid_seq = line[35:-2]
                gene.append(amino_acid_seq)
        elif is_amino_acid_seq:
            if line[-2] != '"':
                amino_acid_seq += line[21:-1]
            else:
                amino_acid_seq += line[21:-2]
                is_amino_acid_seq = False
                gene.append(amino_acid_seq)
    if len(gene) == 3:
         retval.append(gene)
    f.close()
    return retval


if __name__ == "__main__":
    import doctest
    doctest.testmod()
