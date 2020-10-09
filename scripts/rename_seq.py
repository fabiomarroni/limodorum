#!/usr/bin/env python

"""

Author: Fabio Marroni
Contact: marroni@appliedgenomics.org
Date of creation: 2018-12-06
Version: 1.1.0

"""

#PProbably one of the worst python scripts ever... Use at your own risk. Fabio.

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='rename_seq',
                                     description='Given a fasta file that has been classified with Kraken '
                                                 'Add the classification to the name of the sequence (removing the additional trinity statistics) '
                                                 'The script has been tested with a fasta file obtained by trinity and classified with kraken. ' 
                                                 'Different formats of the fasta sequence names might or might not work')
    parser.add_argument('--fasta_in', 
                        dest='fasta_in',
                        required=True,
                        help='Fasta file that has been sent to classification')
    parser.add_argument('--kraken_in', 
                        dest='kraken_input',
                        required=True,
                        help='Output of kraken classification. Use the file showing assignment for each read.')
    parser.add_argument('--fasta_output', 
                        dest='fasta_output',
                        required=True,
                        help='Name of the fasta file with the new names.')
    args = parser.parse_args()
    print(args.fasta_in)


fasta= open(args.fasta_in)
newnames= open(args.kraken_input)
newfasta= open(args.fasta_output, 'w')


fastalines=fasta.readlines()
genelines=newnames.readlines()
for line in fastalines:
    if line.startswith('>'):
        myname=line.split(">")[1]
        myname=myname.split()[0]
        mycount=0
        for newname in genelines:
            genename=newname.split()[1]
            orgname=newname.split("\t")[2:3]
            if(myname==genename):
                orgname=" ".join(orgname)
                fullname=">"+genename+" "+orgname+"\n"
                print "Fastaname=",myname,"Genename=",genename,"Orgname=",orgname
                print fullname
                newfasta.write(fullname)
                mycount=mycount+1
                break
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()


