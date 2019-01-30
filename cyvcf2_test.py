

from cyvcf2 import VCF

for variant in VCF('some.vcf.gz'): # or VCF('some.bcf')
	variant.REF, variant.ALT # e.g. REF='A', ALT=['C', 'T']


	variant.CHROM, variant.start, variant.end, variant.ID, \
				variant.FILTER, variant.QUAL

	# numpy arrays of specific things we pull from the sample fields.
	# gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
	variant.gt_types, variant.gt_ref_depths, variant.gt_alt_depths # numpy arrays
	variant.gt_phases, variant.gt_quals, variant.gt_bases # numpy array


	## INFO Field.
	## extract from the info field by it's name:
	variant.INFO.get('DP') # int
	variant.INFO.get('FS') # float
	variant.INFO.get('AC') # float

	# convert back to a string.
	str(variant)


	## sample info...

	# Get a numpy array of the depth per sample:
    dp = variant.format('DP')
    # or of any other format field:
    sb = variant.format('SB')
    assert sb.shape == (n_samples, 4) # 4-values per

# to do a region-query:

vcf = VCF('some.vcf.gz')
for v in vcf('11:435345-556565'):
    if v.INFO["AF"] > 0.1: continue
    print(str(v))


