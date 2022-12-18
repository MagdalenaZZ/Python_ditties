#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
import sys
import os.path
import argparse
#import csv
import pysam
import re
#from subprocess import call
import subprocess
from pprint import pprint

"""

Script for parsing a small CNV vcf and outputting stats of regions covered by BED

"""


epi = ('\
    \n\
	File parser, allowing counting of CNV variant from VCF files\n\
    \n\
')

# Describe what the script does
parser = argparse.ArgumentParser(description='This script parses a VCF file and extracts regions from a BED file', epilog= epi, formatter_class=argparse.RawTextHelpFormatter)

# Get inputs
parser.add_argument('-i', '--input', default=None, dest='vcf', nargs='+',action='append', type=str, help="VCF files")
parser.add_argument('-b', '--bed', default=None, dest='bed', action='store', required=True, help="BED file")

# Check for no input
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args  = parser.parse_args()

# Check if input files exist
if not os.path.isfile(args.bed)==True:
    print("Cannot find input file ",args.bed)
    sys.exit(1)


#print ("BED file: ",args.bed)
#print ("VCF files: ", args.vcf)


# Create summary output

#print (args.vcf[0][0])
vcf=args.vcf[0][0]
output=''.join([args.bed,'.ovls'])
ou=open(output, 'w')

# Create match output
output=''.join([args.bed,'.rp'])
rp=open(output, 'w')


# Parse the input BED
bed= open(args.bed, 'r')

# Get the total length of reference
bedtot=0

for entry in bed:
    ent = entry.split("\t")
    #print (ent[2],ent[1],int(ent[2])-int(ent[1])+1,sep='\t')
    bedtot+=(int(ent[2])-int(ent[1])+1)

bed.close()


## For each VCF file, do bedtools intersect

# Make output dict
o={}

###############
def print_dict(dct):
    for item, amount in dct.iteritems():
        print('{}:{}'.format(item,amount),',',end='',sep='',file=ou)
    print('\t',end='',file=ou)
################

# Header
print("File","Feature","NoFragsCov","FragsCov","Events","NoEventsPASS","EventsPASS", \
      "EventPass%","TotLen","TotLenPASS", "TotLenPASS%","CovMinPass","CovMaxPass","CovMinAny","CovMaxAny","CovAveALL","CovAvePASS","SharedFilt%","SharedEvent%","SharedCov%","AveDiffCov","AveCov",sep='\t',file=ou)




# Get comparison values for everything

args.vcf=args.vcf[0]

for vcf in args.vcf:
    #print (vcf)
    #["bedtools", "intersect","-b",args.bed,"-a",vcf,"-wo"])
    proc = subprocess.Popen(["bedtools", "intersect","-b",args.bed,"-a",vcf,"-wo"], stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    #print ("OUT:", out)
    out = out.split("\n")

    # For each line in the output, pick out those that overlap
    # Save the BED range for each exon
    for el in out:
        # remove empty lines
        if not (re.match('\w+', el)):
            #print ("Match: ", el,":")
            next
        # Process working lines
        else:
            ele=el.split("\t")
            #print(ele)
            maxstart=max(int(ele[1]),int(ele[13]))
            minend=min(int(ele[2]),int(ele[14]))
            flen=int(minend)-int(maxstart)+1
            #print(maxstart,ele[1],ele[13],minend,ele[2],ele[14],flen,sep='\t')
            bedfeature=''.join([ele[12],"_",ele[13],"_",ele[14]])
            bedregion=''.join([ele[0],"_",ele[1],"_",ele[2]])
            event=ele[3].split(":")[1]
            eventR=ele[15].split(":")[1]

            if vcf not in o:
                o[vcf] = {}
            if ele[23] not in o[vcf]:
                o[vcf][ele[23]]={}
            if bedfeature not in o[vcf][ele[23]]:
                o[vcf][ele[23]][bedfeature]={}
            if bedregion not in o[vcf][ele[23]][bedfeature]:
                o[vcf][ele[23]][bedfeature][bedregion]={}
            o[vcf][ele[23]][bedfeature][bedregion]['Len']=flen
            o[vcf][ele[23]][bedfeature][bedregion]['Event'] =event
            o[vcf][ele[23]][bedfeature][bedregion]['Filter'] =ele[6]
            o[vcf][ele[23]][bedfeature][bedregion]['Copies'] =ele[9]
            o[vcf][ele[23]][bedfeature][bedregion]['Name'] =ele[11]
            o[vcf][ele[23]][bedfeature][bedregion]['EventR'] = eventR
            o[vcf][ele[23]][bedfeature][bedregion]['FilterR'] = ele[18]
            o[vcf][ele[23]][bedfeature][bedregion]['CopiesR'] =ele[21]

            #print("Events",ele[23],bedfeature,bedregion)




# Summarise all PASS canvas fragments for the same feature
# Summarise all canvas fragments for the same feature
# Max depth PASS




def calculate_fragments(data):
    ''' Calculate how many fragments in CNV caller covers REF feature '''
    frags={}

    # Summarise and calculate for each "exon" in "gene"
    for bedfeature, value2 in o[vcf][feature].items():

        # How many canvas fragments cover the same feature?
        for key in o[vcf][feature][bedfeature].keys():
            if key in frags:
                frags[key]+=1
            else:
                frags[key]=1
        # Output all canvas fragments for the same feature
        cnvres = o[vcf][feature][bedfeature].keys()
        cnvres = ','.join(cnvres)
        # print(cnvres,'\t',end='',sep='\t',file=ou)

    print(len(frags),sep='\t',end='\t',file=ou) #3
    print_dict(frags) #4


def calculate_events(data):

    '''Tally up gain, loss, and PASS/fail across all regions of the feature'''

    # How large % of the feature is covered by PASS gain/loss/loh,ref?
    totlen=0
    totlenpass=0
    events={}
    eventspass={}
    eventlen={}


    # Are there both gain and loss in the same feature?

    for bedfeature, value2 in data.items():

        for bedregion, value3 in data[bedfeature].items():
            cops = int(data[bedfeature][bedregion]['Copies'])
            # print("1", bedregion, cops)

            # if data[bedfeature][bedregion]['Event'] == 'Gain':
            if data[bedfeature][bedregion]['Event'] in events:
                events[data[bedfeature][bedregion]['Event']] += 1
            else:
                events[data[bedfeature][bedregion]['Event']] = 1
                #print("New event:",data[bedfeature][bedregion]['Event'])

            if data[bedfeature][bedregion]['Filter'] == 'PASS':
                if data[bedfeature][bedregion]['Event'] in eventspass:
                    eventspass[data[bedfeature][bedregion]['Event']] += 1
                else:
                    eventspass[data[bedfeature][bedregion]['Event']] = 1

    #print(len(events),end='\t',sep='\t',file=ou)
    print_dict(events) #5
    print(len(eventspass), end='\t', sep='\t',file=ou) #6
    print_dict(eventspass) #7

    # How large part of the total length is covered by gain, loss, etc
    for bedfeature, value2 in data.items():

        for bedregion,value3 in data[bedfeature].items():
            cops = int(data[bedfeature][bedregion]['Copies'])
            lens = int(data[bedfeature][bedregion]['Len'])
            evs = data[bedfeature][bedregion]['Event']
            totlen = totlen + lens
            #print("2", bedregion, lens, cops)

            # Only count PASS events
            if data[bedfeature][bedregion]['Filter'] == 'PASS':
                totlenpass = totlenpass + lens

                if evs in eventlen:
                    eventlen[evs] += lens
                    #print ("Is",evs,eventlen[evs])
                else:
                    eventlen[evs]=lens
                    #print("St", evs, eventlen[evs])

    for even,eventlens in eventlen.items():
        #print(even, totlen, (eventlen[even]/totlen))
        eventlen[even]=(eventlen[even]/totlen)

    #print("SEP",totlenpass/totlen,sep='\t',end='\t',file=ou)
    #print_dict(eventlen)
    print_dict(eventlen) #8
    print(totlen, totlenpass,end='\t',sep='\t',file=ou) #9
    return(totlen,totlenpass) #10


def calculate_coverage(data,totlen):

    '''Calculate minimum, maximum and average coverage across all regions of the feature'''

    # Max min and average coverage
    minpass = int(1000000)
    maxpass = 0
    minany = int(1000000)
    maxany = 0
    dpave = {}
    dpavepass = {}
    dpav = 0
    dpavpass = 0

    # Check max,min and average coverage for the entire feature

    for bedfeature, value2 in data.items():

        for bedregion,value3 in data[bedfeature].items():
            cops = int(data[bedfeature][bedregion]['Copies'])
            lens = int(data[bedfeature][bedregion]['Len'])
            #print("3",bedregion,lens, cops)

            if data[bedfeature][bedregion]['Filter']=='PASS':
                if cops > maxpass:
                    maxpass=cops
                if cops < minpass:
                    minpass=cops

                # Add value to average
                if cops in dpavepass:
                    dpavepass[cops] += lens
                    #print ("Is",cops,dpavepass[cops])
                else:
                    dpavepass[cops]=lens
                    #print("St",cops,dpavepass[cops])
            #else:
                #print ("ANY",data[bedfeature][bedregion]['Filter'])

            # Update man and mean
            if  cops > maxany:
                maxany=cops
            if cops < minany:
                minany = cops

            # Add value to average
            if cops in dpave:
                dpave[cops] += lens
                #print ("Is",cops,dpave[cops])
            else:
                dpave[cops]=lens
                #print("St",cops,dpave[cops])

            # Calculate average for all
            for cove,val4 in dpave.items():
                dpav+=(cove*val4)
            dpav=dpav/totlen

            # Calculate average for PASS
            for cove,val4 in dpavepass.items():
                dpavpass+=(cove*val4)
            if totlenpass > 0:
                dpavpass=dpavpass/totlenpass
            else:
                dpavpass=0

    print((totlenpass / totlen), sep='\t', end='\t', file=ou) # 11

    if minany==int(1000000):
        minany='NA'
    if minpass==int(1000000):
        minpass='NA'

    # Max depth all
    # Average depth PASS
    # Average depth all
    print(minpass, maxpass, minany, maxany, "%0.2f" % dpav, "%0.2f" % dpavpass, sep='\t', end='\t', file=ou) # 12,13,14,15,16,17



def calculate_recall(data):

    ''' This section is for comparing reference file -b with inputs -i, if there is match in event, length and coverage'''


    # Check how many % of any region is shared, and % of PASS region shared event
    #print("DATA:",data)

    for bedfeature, value2 in data.items():

        shared_filt=0
        shared_event=0
        shared_cov=0
        totlen = 0
        totlenpass = 0
        diff=0
        covs=0

        for bedregion,value3 in data[bedfeature].items():

            #try:
            event = data[bedfeature][bedregion]['Event']
            #except:
            #    print("EXCEPT:",vcf, ele[23], bedfeature, bedregion)
            #try:
            eventR = data[bedfeature][bedregion]['EventR']
            #except:
                #print("EXCEPT:",vcf, ele[23], bedfeature, bedregion)

            filt = data[bedfeature][bedregion]['Filter']
            filtR = data[bedfeature][bedregion]['FilterR']

            cops = int(data[bedfeature][bedregion]['Copies'])
            copsR = int(data[bedfeature][bedregion]['CopiesR'])

            lens = int(data[bedfeature][bedregion]['Len'])

            # Calcualate lengths
            totlen = totlen + lens
            #print(lens)
            # Only count PASS events
            if filt == 'PASS':
                totlenpass = totlenpass + lens



            # Are they shared PASS/FAIL %
            if filt==filtR:
                shared_filt+=lens
                #print (shared_filt, filt, filtR)
            #else:
                #print("ELSE",shared_filt, filt, filtR)
                #pass

            # Are they shared gain/loss %
            if event==eventR:
                shared_event += lens

            # With a secret save for LOH
            if (event=="LOH" and eventR=="LOSS") or (event=="LOSS" and eventR=="LOH"):
                shared_event += lens

            # Are they shared copy-number %
            if cops==copsR:
                shared_cov += lens
                #print ("IS", bedregion,lens,cops, copsR,sep='\t')
            #else:
                #print("NOT",bedregion,lens,cops, copsR, sep='\t')

            # Average coverage difference
            diff+=((abs(int(cops)-int(copsR)))*lens)
            #print((abs(cops-copsR)), cops, copsR, cops-copsR ,lens,sep='\t')
            #if data[bedfeature][bedregion]['Filter'] == 'PASS':

            # Average coverage
            covs+=(int(cops)*lens)



    # Calculate for the feature, if there is match
        print("%0.2f" % (shared_filt/totlen), "%0.2f" % (shared_event/totlen), "%0.2f" % (shared_cov/totlen), "%0.2f" % (diff/totlen), "%0.2f" %  (covs/totlen) , sep='\t', end='\t',file=ou) # 18,19,20,21,22

        return(totlen,(shared_filt),(shared_event),(shared_cov),(diff))


#####################################################################
#### Outputs ####

print("File","TotLen","FiltOvl%","EventOvl%","CovOvl%","CovDiff%",sep='\t',file=rp)

for vcf,value in o.items():

    print(vcf,'\t',sep='\t',end='',file=rp)

    totlenSum=0.0
    filtSum=0.0
    eventSum=0.0
    covSum=0.0
    diffSum=0.0

    for feature, value in o[vcf].items():
        print(vcf,feature, sep='\t', end='\t', file=ou)  # 1,2

        calculate_fragments(o[vcf][feature])
        #print ("LEN1:",o[vcf][feature])
        totlen, totlenpass = calculate_events(o[vcf][feature])
        #print ("LEN2:",o[vcf][feature])
        calculate_coverage(o[vcf][feature],totlen)
        #print ("LEN3:",o[vcf][feature])
        totlen,shared_filt,shared_event,shared_cov,diffs= calculate_recall(o[vcf][feature])
        #print ("LEN4:",o[vcf][feature])

        totlenSum+=totlen
        filtSum+=shared_filt
        eventSum+=shared_event
        covSum+=shared_cov
        diffSum+=diffs

        print('\n', end='', file=ou)

    print ( bedtot, totlenSum, "%0.4f" % (totlenSum/bedtot), "%0.4f" % (filtSum/bedtot), "%0.4f" % (eventSum/bedtot),  "%0.4f" % (covSum/bedtot),  "%0.4f" % (diffSum/bedtot), sep='\t', end='\n', file=rp)
    #print (  "%0.4f" %  (totlenSum/bedtot), sep='\t', end='\n', file=rp)

ou.close()
rp.close()


quit()





