from pysam import VariantFile


bcf_in = VariantFile("test.vcf.gz")  # auto-detect input format
bcf_out = VariantFile('test.vcf.pysam.vcf.gz', 'w', header=bcf_in.header)


# Add the HP field to header. Say its a string and can take any no. of values. It depends what format you want to give.
bcf_in.header.formats.add("HP",".","String","Some description")


for rec in bcf_in.fetch():
    print (rec.format.keys())
    #for sample in rec.samples:
    rec.samples[0]['HP'] = "Hola"

    #bcf_out.write(rec)




"""

#read the input file
myvcf=VariantFile("/Users/magz/test.vcf.gz","r")
vcf_out = VariantFile("/Users/magz/test.pysam.vcf.gz", 'wu', header=myvcf.header)


# Add the HP field to header. Say its a string and can take any no. of values. It depends what format you want to give.
myvcf.header.formats.add("HP",".","String","Some description")


# An example showing how to add 'HP' information. Please explore further.
for variant in myvcf.fetch():
    print (variant.info)
    #for sample in variant.samples:
        #variant.samples[sample]['HP']="Hola"
    #print variant

    #vcf_out.write(variant)


exit(0)


"""




