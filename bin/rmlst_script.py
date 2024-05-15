#!/usr/bin/env python

# Upload contigs file to PubMLST rMLST species identifier via RESTful API
# Written by Keith Jolley
# Copyright (c) 2018, University of Oxford
# Licence: GPL3

import sys, requests, argparse, base64
parser = argparse.ArgumentParser()
parser.add_argument('--files', '-f', type=str, nargs='+', required=True, help='assemblies (FASTA format)')
args = parser.parse_args()

def main():
    uri = 'http://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1/sequence'
    fasta_input = args.files
    print("Sample_name\tRank\tTaxon\tSupport\tTaxonomy")
    for sequence in fasta_input:
        with open(sequence, 'r') as x: 
            fasta = x.read()
        payload = '{"base64":true,"details":true,"sequence":"' + base64.b64encode(fasta.encode()).decode() + '"}'
        response = requests.post(uri, data=payload)
        if response.status_code == requests.codes.ok:
            data = response.json()
            try:
                data['taxon_prediction']
                for match in data['taxon_prediction']:
                    print(sequence + "\t" + match['rank'] + "\t" + match['taxon'] + "\t" + str(match['support']) + "%" + "\t" + match['taxonomy'])
                    # print("Rank: " + match['rank'])
                    # print("Taxon:" + match['taxon'])
                    # print("Support:" + str(match['support']) + "%")
                    # print("Taxonomy" + match['taxonomy'] + "\n")
            except KeyError:
                print("No match for " + sequence)
                
        else:
            print(response.text)

if __name__ == "__main__":
    main()
