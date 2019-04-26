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
