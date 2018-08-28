import sys
import multiprocessing
import pybedtools

# get example GFF and BAM filenames
gff = pybedtools.example_filename('gdc.gff')
bam = pybedtools.example_filename('gdc.bam')


# Some GFF files have invalid entries -- like chromosomes with negative coords
# or features of length = 0.  This line removes them and saves the result in a
# tempfile
g = pybedtools.BedTool(gff).remove_invalid().saveas()

# This function will eventually be run in parallel, applying the filter above
# to several different BedTools simultaneously
def subset_featuretypes(featuretype):
    result = g.filter(featuretype_filter, featuretype).saveas()
    return pybedtools.BedTool(result.fn)

# This function performs the intersection of a BAM file with a GFF file and
# returns the total number of hits.  It will eventually be run in parallel.
def count_reads_in_features(features_fn):
    """
    Callback function to count reads in features
    """
    # BAM files are auto-detected; no need for an `abam` argument.  Here we
    # construct a new BedTool out of the BAM file and intersect it with the
    # features filename.

    # We use stream=True so that no intermediate tempfile is
    # created, and bed=True so that the .count() method can iterate through the
    # resulting streamed BedTool.
    return pybedtools.BedTool(bam).intersect(
                             b=features_fn,
                             stream=True).count()


# Set up a pool of workers for parallel processing
pool = multiprocessing.Pool()

# Create separate files for introns and exons, using the function we defined
# above
featuretypes = ('intron', 'exon')
introns, exons = pool.map(subset_featuretypes, featuretypes)

for int in introns:
    print(int)



