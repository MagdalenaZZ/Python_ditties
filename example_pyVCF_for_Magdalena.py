import vcf

# reading the VCF file
vcf_reader = vcf.Reader(filename=vcf_input)

    for i, r in enumerate(vcf_reader):

        try:

        	# hash_fields contains all the information from INFO fields
            hash_fields = dict(r.INFO)

            # r.samples contains the name or names included in the VCF file
            # r.FORMAT contains the FORMAT info

            # here I update the FORMAT info together with the INFO info
            hash_fields.update(dict(zip(r.samples[0].data._fields, r.samples[0].data)))




        except:

            raise RuntimeError('<NAME_OF_FUNCTION>: Some error has occurred in variant line')




# First you define the header
target_reader.infos['SR'] = _Info('SR',0,'Float',"Percentage of the samples within the run with the variant")

# Then you enrich that field what whatever you want
r_target.INFO['SR'] = 5


r.FORMAT['NEW'][0]

