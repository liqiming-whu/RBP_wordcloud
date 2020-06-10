#!/usr/bin/env python3
import os
from Bio import Entrez, Medline
from wordcloud import WordCloud
import matplotlib.pyplot as plt


Entrez.email = "liqiming1914658215@gmail.com"
Entrez.api_key = "c80ce212c7179f0bbfbd88495a91dd356708"

stopwords = set(line.rstrip() for line in open("stopwords.txt"))
add = set(line.rstrip() for line in open("stopwords_add.txt"))
stopwords = stopwords | add


def get_count(database, term):
    handle = Entrez.egquery(term=term)
    record = Entrez.read(handle)
    for row in record["eGQueryResult"]:
        if row["DbName"] == database:
            count = row["Count"]
    return count


def search(database, keywords, count):
    handle = Entrez.esearch(db=database, term=keywords, retmax=count)
    record = Entrez.read(handle)
    return record["Count"], record["IdList"]


def get_abstract(database, idlist):
    text = ""
    handle = Entrez.efetch(db=database, id=idlist, rettype="medline", retmode="text")
    records = Medline.parse(handle)
    for record in records:
        try:
            pmid = str(record.get("PMID", "?"))
            print("pmid:", pmid)
            abstract = record.get("AB", "?")
            text += abstract
        except:
            continue
    return text


def save_text(text):
    with open("result.txt", "w", encoding="utf-8") as f:
        f.write(text)


def wordcloud(text):
    """
    params:
    stopwords:set, default:None
    collections:Include binary phrases or not
    """
    wordcloud = WordCloud(background_color="white", stopwords=stopwords,
                          collocations=False, width=1000, height=750, margin=2).generate(text)
    plt.imshow(wordcloud)
    plt.axis("off")
    plt.show()
    wordcloud.to_file("wordcloud.png")


def main():
    query = "bladder cancer OR bladder tumor OR bladder carcinoma OR carcinoma of bladder OR urothelial carcinoma"
    if(os.path.exists("result.txt")):
        text = open("result.txt").read()
        print("results laoded.")
    else:
        count = get_count("pubmed", query)
        print("count:", count)
        _, idlist = search("pubmed", query, count)
        text = get_abstract("pubmed", idlist)
        save_text(text)
    wordcloud(text)


if __name__ == '__main__':
    main()


