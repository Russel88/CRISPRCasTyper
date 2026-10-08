import logging
import sys

import pandas as pd
import pyrodigal_gv


class Prodigal(object):

    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run_prod(self) -> None:

        if not self.redo:
            logging.info("Predicting ORFs with pyrodigal-gv")

            meta_mode = self.prod == "meta"
            vf = pyrodigal_gv.ViralGeneFinder(meta=meta_mode)
            trained = meta_mode

            proteins_handle = open(self.prot_path, "w")
            genes_table = []

            # In single-genome mode pyrodigal-gv requires training
            #  before calling find_genes.
            if not trained:
                try:
                    first_seq = next(iter(self.seq_dict.values()))
                except StopIteration:
                    logging.critical(
                        "No sequences available for gene prediction"
                    )
                    sys.exit()
                vf.train(bytes(str(first_seq), encoding="utf-8"))
                trained = True

            for contig_idx, (contig, seq) in enumerate(
                self.seq_dict.items(), start=1
            ):
                genes = vf.find_genes(bytes(str(seq), encoding="utf-8"))
                if not genes:
                    logging.warning("No genes predicted for contig %s", contig)
                    continue

                for gene_idx, gene in enumerate(genes, start=1):
                    header = f">{contig}_{gene_idx} # {gene.begin} #"
                    f" {gene.end} # {gene.strand}"
                    prot = gene.translate(include_stop=False)
                    proteins_handle.write(header + "\n")
                    proteins_handle.write(prot + "\n")
                    genes_table.append(
                        (contig, gene.begin, gene.end, gene.strand, gene_idx)
                    )

            proteins_handle.close()

            if not genes_table:
                logging.critical("Pyrodigal-gv failed to predict any genes")
                sys.exit()

            # Make gene table
            self.genes = pd.DataFrame(
                genes_table,
                columns=("Contig", "Start", "End", "Strand", "Pos"),
            )
            self.genes.to_csv(self.out + "genes.tab", index=False, sep="\t")

    def get_genes(self) -> None:
        # Genes table is constructed during pyrodigal-gv prediction
        if not hasattr(self, "genes"):
            try:
                self.genes = pd.read_csv(self.out + "genes.tab", sep="\t")
            except Exception:
                logging.error(
                    "Gene predictions not found; run pyrodigal-gv first."
                )
                sys.exit()
