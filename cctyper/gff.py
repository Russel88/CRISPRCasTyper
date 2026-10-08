from __future__ import annotations

import logging
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


class GFF(object):

    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def get_genes(self) -> None:

        def unwrap_attributes(df: pd.DataFrame) -> pd.DataFrame:
            attribute_names = set()
            for attributes in df["attributes"]:
                attributes_list = [n for n in attributes.split(";") if n]
                for attr in attributes_list:
                    attribute_name = attr.split("=")[0]
                    attribute_names.add(attribute_name)

            for attribute_name in attribute_names:
                df = df.assign(**{attribute_name: None})

            def extract_attribute_values(row):
                attributes_list = [
                    n for n in row["attributes"].split(";") if n
                ]
                for attr in attributes_list:
                    attribute_name, attribute_value = attr.split("=")
                    row[attribute_name] = attribute_value
                return row

            df = df.apply(lambda row: extract_attribute_values(row), axis=1)
            df = df.drop(columns=["attributes"])
            return df

        custom_header = [
            "Contig",
            "source",
            "type",
            "Start",
            "End",
            "score",
            "Strand",
            "phase",
            "attributes",
        ]
        gff = pd.read_csv(
            self.gff, sep="\t", comment="#", names=custom_header, engine="c"
        )
        gff = gff[gff["type"] == "CDS"]
        if gff.empty:
            logging.error("No CDS features found in GFF %s", self.gff)
            sys.exit()

        gff["Strand"] = (
            gff["Strand"].map({"+": 1, "-": -1}).fillna(0).astype(int)
        )
        gff["Start"] = pd.to_numeric(gff["Start"], errors="coerce")
        gff["End"] = pd.to_numeric(gff["End"], errors="coerce")
        gff = gff.dropna(subset=["Start", "End"])
        if gff.empty:
            logging.error(
                "GFF file missing valid Start/End coordinates after filtering"
            )
            sys.exit()
        gff[["Start", "End"]] = gff[["Start", "End"]].astype(int)

        gff = unwrap_attributes(gff)
        prot_records = self._load_proteins()
        prot_ids = set(prot_records.keys())

        protein_id = None
        if "protein_id" in gff.columns:
            protein_id = gff["protein_id"]
        elif "ID" in gff.columns:
            ids = gff["ID"]
            if ids.str.startswith("cds-").all():
                protein_id = ids.str.split("-", 1).str[1]
            elif ids.str.match(r"^\D+_\d+$").all():
                protein_id = ids
            elif ids.str.match(r"^\d+_\d+$").all():
                protein_id = gff["Contig"] + "_" + ids.str.split("_").str[1]
            elif set(ids).issubset(prot_ids):
                protein_id = ids
            else:
                logging.error(
                    "GFF file does not match expected format. Please check"
                    " that it meets the format requirements."
                )
                sys.exit()
        else:
            logging.error("No protein_id or ID column found in GFF.")
            sys.exit()

        gff["protein_id"] = protein_id
        gff = gff[["Contig", "Start", "End", "Strand", "protein_id"]].dropna(
            subset=["protein_id"]
        )
        gff = gff.sort_values(["Contig", "Start"])
        gff["Pos"] = gff.groupby("Contig").cumcount() + 1

        genes = gff[["Contig", "Start", "End", "Strand", "Pos", "protein_id"]]
        genes.to_csv(self.out + "genes.tab", index=False, sep="\t")

        out_faa = Path(self.prot_path)
        with out_faa.open("w") as fh:
            for _, row in genes.iterrows():
                pid = str(row["protein_id"])
                rec = prot_records.get(pid)
                if rec is None:
                    logging.error(
                        "Protein %s referenced in GFF not found in"
                        " protein FASTA",
                        pid,
                    )
                    sys.exit()
                header = f">{pid} # {row['Contig']}_{int(row['Pos'])} #"
                f" {int(row['Start'])} # {int(row['End'])} #"
                f" {int(row['Strand'])}"
                fh.write(f"{header}\n{str(rec.seq)}\n")

    def _load_proteins(self):
        records = {}
        for record in SeqIO.parse(self.prot, "fasta"):
            records[record.id] = record
        if not records:
            logging.error("No protein sequences found in %s", self.prot)
            sys.exit()
        return records
