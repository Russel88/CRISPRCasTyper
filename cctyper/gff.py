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
                attributes_list = [n for n in attributes.split(';') if n and n!='']
                for attr in attributes_list:
                    attribute_name = attr.split('=')[0]
                    attribute_names.add(attribute_name)

            # Create new columns with attribute names and default values
            for attribute_name in attribute_names:
                df = df.assign(**{attribute_name: None})

            # Function to extract attribute values and assign them to columns
            def extract_attribute_values(row):
                attributes_list = [n for n in row['attributes'].split(';') if n and n!=''] 
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

        gff = pd.read_csv(self.gff,sep="\t",comment="#",names=custom_header,engine="c")
        gff = gff[gff['type'] == 'CDS']
        gff['Strand'] = gff['Strand'].map({'+': 1, '-': -1})

        gff=unwrap_attributes(gff)
        cols = [col for col in gff.columns if col in [n for n in custom_header if n[0].isupper()]]
        if "ID" in gff.columns:
            if gff['ID'].str.startswith('cds-').all():
                gff["protein_id"]=gff["ID"].str.split("-").apply(lambda x: x[1])
            elif gff["ID"].str.match(r'^\D+_\d+$').all():
                gff["protein_id"]=gff["ID"]
            elif gff["ID"].str.match(r'^\d+_\d+$').all():
                gff["protein_id"] = gff["Contig"] + "_" +gff["ID"].str.split("_").apply(lambda x: x[1])
            else:
                prot_ids = subprocess.run(['sed', '-n', 's/^>//p', self.prot_path],capture_output=True).stdout.decode().splitlines()
                if set(prot_ids).issubset(gff['ID']):
                    gff["protein_id"]=gff["ID"]
                else:
                    logging.error("GFF file does not match expected format. Please check that it meets the format requirements.")
                    sys.exit()
        gff=gff[cols+["protein_id"]].dropna(subset=["protein_id"])
        gff = gff.sort_values(['Contig', 'Start'])
        gff['Pos'] = gff.groupby('Contig').cumcount() + 1
        self.genes = gff[cols+["Pos","protein_id"]]
        self.genes.to_csv(self.out+'genes.tab', index=False, sep='\t')