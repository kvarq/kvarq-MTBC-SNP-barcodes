
VERSION = '0.1'
GENES_COMPATIBILITY = '0.0'

import os.path

from kvarq.genes import Genome, Reference, SNP, Test, Testsuite, Genotype

def tsv2SNPs(path, genome, reference):

    tests = []
    for line in file(path):

        parts = line.strip().split('\t')
        name = parts[0]
        pos = int(parts[1])
        bases = parts[2].split('/')

        snp = SNP(genome=genome, pos=pos, orig=bases[0], base=bases[1])
        test = Test(snp, Genotype(name), reference)
        tests.append(test)

    return tests

here = os.path.dirname(__file__)
genome = Genome(os.path.join(here, 'MTB_ancestor_reference.bases'), 'MTB ancestor')
coll14 = Reference('Coll et al (2014) -- doi: 10.1038/ncomms5812')
SNPs = tsv2SNPs(os.path.join(here, 'coll14.tsv'), genome, coll14)

coll14 = Testsuite(SNPs, VERSION)

