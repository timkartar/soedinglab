#/usr/bin/env python

''' Parse dosage data from VCF file.
    Usage:
       from iotools.readvcf import ReadVCF
       readvcf = ReadVCF(vcf_filepath)
       snpinfo = readvcf.snpinfo
          ...
    Returns:
       a) information of all the variants read
            - chrom
            - pos
            - id
            - reference allele
            - alternate allele
            - minor allele frequency
       b) genotype matrix
       c) identity of the donors
'''

import numpy as np
import gzip

import collections

SNPINFO_FIELDS = ['chrom', 'varid', 'bp_pos', 'ref_allele', 'alt_allele', 'maf']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()

SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

GENEINFO_FIELDS = ['name', 'ensembl_id', 'chrom', 'start', 'end']
class GeneInfo(collections.namedtuple('_GeneInfo', GENEINFO_FIELDS)):
    __slots__ = ()

class ReadVCF:


    _read_genotype_once = False


    def __init__(self, filepath, mode="DS"):
        self._filepath = filepath
        self._mode = mode


    @property
    def snpinfo(self):
        self._run_once()
        return tuple(self._snpinfo)


    @property
    def dosage(self):
        self._run_once()
        return self._dosage


    @property
    def donor_ids(self):
        self._run_once()
        return tuple(self._donor_ids)


    def _run_once(self):
        if self._read_genotype_once:
            return
        self._read_genotype_once = True
        self._read_dosage()


    def _read_dosage(self):
        dosage = list()
        snpinfo = list()
        with gzip.open(self._filepath, 'r') as vcf:
            for line in vcf:
                linestrip = line.decode().strip()
                if linestrip[:2] == '##': continue
                if linestrip[:6] == '#CHROM':
                    linesplit = linestrip.split("\t")
                    donor_ids = linesplit[9:]
                else:
                    linesplit = linestrip.split("\t")
                    if len(linesplit[3]) ==1 and len(linesplit[4]) ==1:
                    	if linesplit[4] != SNP_COMPLEMENT[linesplit[3]]:
		                    chrom = int(linesplit[0])
		                    pos   = int(linesplit[1])
		                    varid = linesplit[2]
		                    ref   = linesplit[3]
		                    alt   = linesplit[4]
		                    if self._mode == "DS":
		                        dsindx = linesplit[8].split(':').index("DS")
		                        ds = [x.split(':')[dsindx] for x in linesplit[9:]]
		                        gtindx = linesplit[8].split(':').index("GT")
		                        for i, x in enumerate(ds):
		                            if x == ".":
		                                gt = linesplit[9+i].split(':')[gtindx]
		                                if len(gt) == 3 and gt[0] != "." and gt[2] != ".":
		                                    ds[i] = float(int(gt[0]) + int(gt[2]))

		                    elif self._mode == "GT":
		                        gtindx = linesplit[8].split(':').index("GT")
		                        gt = [x.split(':')[gtindx] for x in linesplit[9:]]
		                        ds = [ float(int(x[0]) + int(x[2])) if len(x) == 3 and x[0] != "." and x[2] != "." else "." for x in gt ]

		                    ds_notna = [float(x) for x in ds if x != "."]
		                    freq = sum(ds_notna) / 2 / len(ds_notna)
		                    if freq <=0.5:
		                    	maf=freq
		                    else:
		                    	maf = 1 - freq
		                    #maf = freq
		                    if maf >= 0.1:#and freq <= 0.9:
			                    snpdosage = [float(x) if x != '.' else 2 * freq for x in ds]
			                    #if freq > 0.5:
			                    #    maf = 1 - freq
			                    #    ref = linesplit[4]
			                    #    alt = linesplit[3]
			                    #    snpdosage = [2 - x for x in snpdosage]
			                    
			                    this_snp = SnpInfo(chrom      = chrom,
			                                       bp_pos     = pos,
			                                       varid      = varid,
			                                       ref_allele = ref,
			                                       alt_allele = alt,
			                                       maf        = maf)

			                    dosage.append(snpdosage)
			                    snpinfo.append(this_snp)

        self._dosage = np.array(dosage)
        self._snpinfo = snpinfo
        self._donor_ids = donor_ids
