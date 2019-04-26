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


def getdata(numberofrecs):
    print('Checking Archives...')
    try:
        readrecordsdf = pd.read_csv(('selected_records' + str(numberofrecs) + '.csv'), index_col=0)
        queryinput = readrecordsdf['Searched'][0]
        print('{} records retrieved for \'{}\''.format(numberofrecs, queryinput))
    except OSError:
        if int(numberofrecs) < 10001:
            [queryinput, idlist] = searchpubmed(return_max=numberofrecs, database="pubmed")
            readrecordsdf = fetchrecord(idlist, numberofrecs, queryinput)
        else:
            print('please check size and format restrictions')
            sys.exit()
    return {'SearchQuery': queryinput, 'RecordsDataFrame': readrecordsdf}

def topicnetwork(inputarray, numtops):
    plt.figure()
    topiclist = []
    for i in range(int(numtops)):
        topiclist.append('Topic ' + str(i))
    G = nx.Graph()
    G.add_nodes_from(topiclist)
    for i in range(int(numtops)):
        getstarseries = []
        getstarseries.append(topiclist[i])
        for term in inputarray['TopicTerms'][i]:
            getstarseries.append(term)
        G.add_star(getstarseries)
    nx.draw_networkx(G, node_size=10, edge_color='r', font_color='k', font_weight='bold')
    plt.title('Topic Term Distribution Network')
    plt.show()













if __name__ == "__main__":
    ## Get Data
    numberofrecs = input("How many listings should we download (100 < recommended < 10,000)? \n")
    numsearchtopics = min(max(10, round(int(numberofrecs) / 200, -1)), 30)  # limits for num topics
    print('\nSeek {} topics in {} documents'.format(int(numsearchtopics), int(numberofrecs)))

    getdatadict = getdata(numberofrecs)
    querysearched = getdatadict['SearchQuery']
    getresultdf = getdatadict['RecordsDataFrame']
    cleanedresultdf = cleandata(getresultdf)

    ##Topic Model:
    customstopwords = stop_words.get_stop_words('en', cache=False)
    newstopwords = ['effect', 'increased', 'increase', 'decreased', 'decrease', 'inhibit', 'result', \
                    'role', 'regulate', 'via', 'associated', 'associate', 'new', 'inhibitor', 'antagonist', 'agonist', \
                    'dependent', 'independent']
    for term in newstopwords:
        customstopwords.append(term)
    for term in querysearched.split(", "):
        customstopwords.append(term)

    cleanedresultdf['TokensTitles'] = tokenizefortopicmodel(cleanedresultdf.TI, stopwords=customstopwords)
    cleanedresultdf['TokensAbstracts'] = tokenizefortopicmodel(cleanedresultdf.AB, stopwords=customstopwords)

    topicmodeldict = topicmodel(cleanedresultdf.TokensTitles, cleanedresultdf.TokensAbstracts,
                                numtopics=numsearchtopics, numwords=10)
    term_dictionary = topicmodeldict['Dictionary']
    topic_terms_df = topicmodeldict['TopicTermsDF']
    doc_similarity_matrix = topicmodeldict['SimilarityMatrix']

    cleanedresultdf['DocTopics'] = topicmodeldict[
        'DocTopProbAbs']  # topics per doc, tuples (topic, probability), ranked, for topics > minimum probability
    cleanedresultdf['DocTopicsTop'] = topicmodeldict['DocTopProbAbsHighest']  # top topic per doc

    ## Inspect
    print(topic_terms_df['TopicTerms'])
    topicnetwork(topic_terms_df, numsearchtopics)

    pubyearhistogram(cleanedresultdf, topicplotted=querysearched)

    listoftopicdfs = []
    for t in range(int(numsearchtopics)):
        topicdocdf = cleanedresultdf[cleanedresultdf.DocTopicsTop == t]
        pubyearhistogram(topicdocdf, "Topic " + str(t))
        listoftopicdfs.append(topicdocdf.ReleaseDate.dt.year)

    topicdistributionhistogram(listoftopicdfs, cleanedresultdf)

    print('\nDone')
