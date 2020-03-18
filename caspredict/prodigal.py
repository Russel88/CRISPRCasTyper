import os
import subprocess
import logging
import sys

class Prodigal(object):
    
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run_prod(self):

        if not self.redo:
            logging.info('Predicting ORFs with prodigal')

            # Run prodigal
            with open(self.out+'prodigal.log', 'w') as prodigal_log:
                subprocess.run(['prodigal', 
                                '-i', self.fasta, 
                                '-a', self.out+'proteins.faa', 
                                '-p', self.prod], 
                                stdout=subprocess.DEVNULL, 
                                stderr=prodigal_log)

            # Check if succesful
            self.check_rerun()

    def check_rerun(self):
        # Check prodigal output
        if os.stat(self.prot_path).st_size == 0:
            if self.prod == 'single':
                logging.warning('Prodigal failed. Trying in meta mode')
                self.prod = 'meta'
                self.run_prod()
            else:
                logging.critical('Prodigal failed! Check the log')
                sys.exit()
        
