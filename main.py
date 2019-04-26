""" Created on 4 April 2019,
Modified on 26 April 2019

@Author: Anoushkrit Goel 
"""
!pip install sys
import sys
!pip install re
import re
import numpy as np

!pip install pprint 
from pprint import pprint

# Gensim
!pip install gensim 
import gensim
import gensim.corpora as corpora 
from gensim.utils import simple_preprocess

from gensim.models import CoherenceModel # Coherence Model may not be used

# spacy for Lemmatization
!pip install spacy
!python -m download spacy en
import spacy

#Plotting Tools
!pip install pyldavis
import pyldavis 
import pyldavis.gensim #crucial step 
import matplotlib.pyplot as plt
%matplotlib inline

# Network Visualisation 
!pip install networkx
import networkx as nx

#Quality Analysis and Testing Purposes
!pip install logging 
import logging
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.ERROR)

!pip install warnings
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

def searchpubmed(return_max=20, database=None, searchqueryinput=None):
    print('Search' + database)
    if searchqueryinput:
        searchquery = searchqueryinput
    handle = Entrez.esearch(db=database, term=searchquery, retmax=return_max)  # retmax is 20 by esearch default
    queryresponse = Entrez.read(handle)
    handle.close()
    print("\n" + queryresponse['Count'] + " results available for \'" + str(searchquery) + "\', returning " + str(
        return_max))
    queryids = (searchquery, queryresponse['IdList'])
    return queryids

def fetchrecord(inputids, numberofrecs, queryinput):
    print('Downloading Pubmed Records...')
    handle = Entrez.efetch("pubmed", id=str(inputids), rettype="medline", retmode="text")
    records = Medline.parse(handle)
    records = list(records)

    recordsdf = pd.DataFrame(records)
    print('\nFeatures available: {}'.format(recordsdf.columns.values.tolist()))
    recordsdf['Searched'] = str(queryinput)
    reckeys = ['PMID', 'TI', 'AB', 'DP', 'PHST', 'Searched']
    print("\nReturning {} features : {} ".format((len(reckeys) - 1), reckeys[0:(len(reckeys) - 1)]))
    recordsselectdf = recordsdf[reckeys]

    recordsdf.to_csv('full_records' + str(numberofrecs) + '.csv')
    recordsselectdf.to_csv('selected_records' + str(numberofrecs) + '.csv')

    handle.close()
    return recordsselectdf
