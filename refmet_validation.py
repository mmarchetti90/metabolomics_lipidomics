#!/usr/bin/env python3

"""
Class for converting metabolite/lipid names to RefMet standard
"""

### ---------------------------------------- ###

class refmet_validation:
    
    """
    Class for converting metabolite/lipid names to RefMet standard
    
    Parameters
    ----------
    compounds : list of strings
        Metabolites/lipids to be validated.
    
    Attributes
    ----------
    compounds : list of strings
        Metabolites/lipids to be validated.
    base_url: string
        Base URL for RefMet request.
    validation_table : pd.DataFrame
        Results of the validation.
    """
    
    def __init__(self, compounds=['glucose']):
        
        self.compounds = list(compounds)
        
        #self.base_url = 'https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_min.php'
        self.base_url = 'https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_minID.php'
    
    ### ------------------------------------ ###
    
    def validate(self):
        
        compounds_dict = {'metabolite_name': '\n'.join(self.compounds)}
        
        refmet_response = requests.post(self.base_url, data=compounds_dict)
        
        if refmet_response.status_code == 200:
            
            header, *validation_table = [rr.split('\t') for rr in refmet_response.text.split('\n')]
            
            self.validation_table = pd.DataFrame(validation_table, columns=header)
            
        
        else:
            
            print(f'ERROR: status code {refmet_response.status_code}')
            
            self.validation_table = pd.DataFrame([])

### ------------------MAIN------------------ ###

import pandas as pd
import requests
