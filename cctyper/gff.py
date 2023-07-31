import os
import subprocess
import logging
import sys
import re

import pandas as pd

class GFF(object):
    
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def get_genes(self):

        def unwrap_attributes(df):
            attribute_names = set()
            for attributes in df['attributes']:
                attributes_list = attributes.split(';')
                for attr in attributes_list:
                    attribute_name = attr.split('=')[0]
                    attribute_names.add(attribute_name)

            # Create new columns with attribute names and default values
            for attribute_name in attribute_names:
                df = df.assign(**{attribute_name: None})

            # Function to extract attribute values and assign them to columns
            def extract_attribute_values(row):
                attributes_list = row['attributes'].split(';')
                for attr in attributes_list:
                    attribute_name, attribute_value = attr.split('=')
                    row[attribute_name] = attribute_value
                return row

            # Apply the function to each row
            df = df.apply(lambda row: extract_attribute_values(row), axis=1)

            # Drop the original "attributes" columns
            df = df.drop(columns=['attributes'])
            return df

        custom_header = ['Contig', 'source', 'type','Start','End','score','Strand','phase','attributes']

        gff = pd.read_csv("genomic.gff",sep="\t",comment="#",names=custom_header,engine="c")
        gff = gff[gff['type'] == 'CDS']
        gff['Strand'] = gff['Strand'].map({'+': 1, '-': -1})

        gff=unwrap_attributes(gff)
        cols = [col for col in gff.columns if col in [n for n in custom_header if n[0].isupper()]]
        gff=gff[cols+["protein_id"]].dropna(subset=["protein_id"])
        gff = gff.sort_values(['Contig', 'Start'])
        gff['Pos'] = gff.groupby('Contig').cumcount() + 1
        self.genes = gff[cols+["Pos","protein_id"]]
        self.genes.to_csv(self.out+'genes.tab', index=False, sep='\t')