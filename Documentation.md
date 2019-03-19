



#Biomedical Text Analytics














 
1. Introduction
Biomedical Text Analytics is another field where researchers and professionals have worked upon. This field based on the past researches proves that it has potential in predicting the package in   which a disease and other diseases affect a patient. Not only this, researchers gain a lot of benefit from this type of approach as they get to know the fields which require research and the undiscovered correlations among various things.  
Topic Modeling is a subfield of Text Mining and Analytics which further is a subfield of Natural Language Processing which comes under Deep Learning. Deep Learning, on the other hand, is a subfield of Machine Learning. Deep Learning has seen a lot of growth in the year 2013 and is increasing at a tremendous. Every single day, an overwhelming number of articles, research papers, and knowledge content is published. Not only research but companies are also coming into the market with reinforced technologies. Artificial Intelligence is now considered as the 4th Industrial Revolution.
1.1 Topic Modeling
 In natural language understanding (NLU) tasks, there is a hierarchy of lenses through which we can extract meaning — from words to sentences to paragraphs to documents. At the document level, one of the most useful ways to understand text is by analyzing its topics. The process of learning, recognizing, and extracting these topics across a collection of documents is called topic modeling [1].
1.1.1 LSA (Latent Semantic Analysis)
Latent semantic analysis (LSA) is a technique in Natural Language Processing distributional semantics, of analyzing relationships between a set of documents and the terms they contain by producing a set of concepts related to the documents and terms. LSA assumes that words that are close in meaning will occur in similar pieces of text (the distributed computing). A matrix containing word counts per paragraph (rows represent unique words and columns represent each paragraph) is constructed from a large piece of text and a mathematical technique called Singular Value Decomposition (SVD) is used to reduce the number of rows while preserving the similarity structure among columns. Paragraphs are then compared by taking the cosine of the angle between the two vectors (or the dot product between the normalizations of the two vectors) formed by any two columns. Values close to 1 represent very similar words while values close to 0 represent very dissimilar words [2].
For example:  if you search for "dogs" and you search results show a document which doesn't contain dog but instead has word "canine" occurring very frequently, LSA is at play.
LSA reduces the dimensions by characterizing similar word vectors to “topics” and then further reducing the redundancy of the vectors carrying the same intuitive meaning attached to it.

1.1.2 PLSA (Probabilistic Latent Semantic Analysis)
PLSA, or Probabilistic Latent Semantic Analysis, uses a probabilistic method instead of SVD to tackle the problem. The core idea is to find a probabilistic model with latent topics that can generate the data we observe in our document-term matrix. In particular, we want a model P (D, W) such that for any document d and word w, P(d,w) corresponds to that entry in the document-term matrix [1].
1.1.3 LDA (Latent Dirichlet Allocation)
A topic model takes a collection of unlabeled documents and attempts to find the structure or topics in this collection. Note that topic models often assume that word usage is correlated with the topic occurrence. You could, for example, provide a topic model with a set of news articles and the topic model will divide the documents in several clusters according to word usage.

Topic models are a great way to automatically explore and structure a large set of documents: they group, or cluster documents based on the words that occur in them. As documents on similar topics tend to use a similar sub-vocabulary, the resulting clusters of documents can be interpreted as discussing different “topics”.

Latent Dirichlet Allocation (LDA) is an example of a probabilistic topic model. 
 
Fig 1.1.3
The model representation for LDA is shown below with Probabilistic approach 
Fig 1.1.3(2)

1.1.4 Word2Vec
Word2Vec, as the name suggests, converts the words in a corpus to vectors and is corpus specific and not documents specific.
 
Fig 1.1.4
1.1.5 lda2vec
lda2vec, in short, is the hybrid of both the algorithms, LDA and word2vec and helps us clearly interpret the data and maintains the flexibility [3].



	  
2. Our Services
We intend to provide a service which can automate the process of finding concomitant and co-morbid diseases with the help of analysing millions of research papers. The amount of research papers can be very large and is extracted from Public Databases by open-source Platform, NCBI.
The Public Databases include PubMed, ClinVar, PMC, MeSH etc.
 
3. Project Phases
3.1 Problem Identification
Problem of finding correlation between varied diseases, symptoms, causes is still not very evident but can be visualised by the help of Chord Diagrams. 
3.2 Literature Review
Literature review, including but not limited to the theory behind the algorithms is provided in the section 1.1 (Topic Modelling). 
3.3  
4. Pipeline
 
5. Methodology
5.1 Data Extraction 
There are multiple ways possible to extract the data. Some of them are listed below for reference.
5.1.1 BioPython 
BioPython[4] is a library specifically made to scrap multiple PMIDs ( A PMID is the unique identifier number used in PubMed for each article) without NCBI servers blocking the user from accessing the data as this type of Data Scraping is usually characterised as an attack on the server. Hence, to protect the valuable information servers usually block or limit the access of articles for that user.
Bio.Entrez
Bio.Entrez is the library which enables the user id specific scraping of data and only allowing us to scrap the relevant information.
5.1.2 Scrapy
Scrapy, an open source and collaborative framework for extracting the data you need from websites, is a web crawler which possesses the ability to scrap almost all data from any website. But this type of scraping can be very unstructured and can also be classified under brute force attack if not used properly.
5.1.3  

5.2 Data Pre-processing
5.2.1 Data Cleaning 
One of the most crucial elements to Data Pre-processing is getting rid of the unnecessary data that is present in the extracted dataset.
5.3 
 
6. Output
6.1 Proposed Output 
7. References
[1] Joyce Xu. Topic Modelling with LSA, PLSA, LDA & lda2Vec.           https://medium.com/nanonets/topic-modeling-with-lsa-psla-lda-and-lda2vec-555ff65b0b05
 https://towardsdatascience.com/lda2vec-word-embeddings-in-topic-models-4ee3fc4b2843
[2] Latent Semantic Analysis 
https://en.wikipedia.org/wiki/Latent_semantic_analysis 
[3] Christopher Moody. word2vec, LDA, and introducing a new hybrid algorithm: lda2vec.            https://www.slideshare.net/ChristopherMoody3/word2vec-lda-and-introducing-a-new-hybrid-alg orithm-lda2vec-57135994
[4] BioPython https://biopython.org/wiki/Documentation
 
Appendix
 

 

 

