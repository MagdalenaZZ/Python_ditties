from pyCGA.opencgarestclients import OpenCGAClient
import sys
import os

input_file = sys.argv[1]

if not os.path.exists(input_file):
    print('Specified input file "{}" could not be found'.format(input_file))
    sys.exit()

lp_nums = [line.rstrip() for line in open(input_file)]


"""
# need this configuration file, contents:
{
    "version": "v1",
    "rest": {
        "hosts": [
           "bio-prod-opencgainternal-haproxy-01.gel.zone/opencga"
        ]
    }
}
"""
oc = OpenCGAClient(configuration='ocga_config.json', user='interpretationpipeline', pwd='kx6SHa26')

# iterate over all studies
studies = [1000000034, 1000000038, 1000000032, 1000000024]

# iterate over all LPs
for lp in lp_nums:

    # check for cohorts - assume details not found
    cohort_id = False
    study = False

    # iterate over both studies - should only really find in 38
    for x in studies:

        # not finding any results using this client throws an error
        try:
            cohort_details = oc.cohorts.search(study=x, samples=lp)
            if len(cohort_details) != 0:
                cohort_id = cohort_details[0][-1]['name']
                study = x
                break

        # expected, ignore
        except:
            continue

    # if cohort details were found, print LP, study, cohort
    if cohort_id:
        print('{} : {} : {}'.format(lp, study, cohort_id))
    else:
        print('cohort not found in either study: {}'.format(lp))
