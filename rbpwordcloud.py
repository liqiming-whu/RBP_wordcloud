#!/usr/bin/env python3
import os
import csv
from collections import Counter
from Bio import Entrez, Medline
from wordcloud import WordCloud
import matplotlib.pyplot as plt


Entrez.email = "liqiming1914658215@gmail.com"
Entrez.api_key = "c80ce212c7179f0bbfbd88495a91dd356708"

rbps = list(line.rstrip() for line in open("rpb_gene.txt"))

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

def save_text(text, pmid):
    with open("download/"+pmid+".txt", "w", encoding="utf-8") as f:
        f.write(text)

def get_abstract(database, idlist):
    handle = Entrez.efetch(db=database, id=idlist, rettype="medline", retmode="text")
    records = Medline.parse(handle)
    for record in records:
        try:
            pmid = str(record.get("PMID", "?"))
            print(pmid)
            abstract = record.get("AB", "?")
            save_text(abstract, pmid)
        except:
            continue


def merge_txt(path):
    with open("results.txt", "w", encoding="utf-8") as f:
        for file in os.listdir(path):
            text = open(os.path.join(path, file)).read()
            f.write(text)


def build_word_dict(text):
    with open(text,encoding="utf-8") as f:
        word_dict = Counter(f.read().split())
    return word_dict


def get_csv(word_dict, path):
    count_list = []
    for gene in rbps:
        if gene == 'TARDBP':
            count_list.append(["TDP43", word_dict[gene] + word_dict["TDP43"] + word_dict["TDP-43"]])
        elif gene == 'FMR1':
            count_list.append(["FMRP", word_dict["FMR1"] + word_dict["FMRP"]])
        else:
            count_list.append([gene, word_dict[gene]])
    count_list.sort(reverse=True, key=lambda x: x[1])
    with open(path,"w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene", "count"])
        for count in count_list:
            if count[1] > 0:
                writer.writerow(count)


def wordcloud(word_dict):
    rpb_dict = dict()
    for gene in rbps:
        if gene == 'TARDBP':
            rpb_dict['TDP43'] = word_dict[gene] + word_dict["TDP43"] + word_dict["TDP-43"]
        elif gene == 'FMR1':
            rpb_dict["FMRP"] = word_dict["FMR1"] + word_dict["FMRP"]
        elif word_dict[gene] > 0:
            rpb_dict[gene] = word_dict[gene]         
    wordcloud = WordCloud(background_color="white", collocations=False, scale=4,
                          width=1000, height=750, margin=2).fit_words(rpb_dict)
    plt.imshow(wordcloud)
    plt.axis("off")
    plt.show()
    wordcloud.to_file("wordcloud.pdf")


def main():
    query = "RNA binding protein AND neuron"
    if(os.path.exists("results.txt")):
        print("results laoded.")
    else:
        count = get_count("pubmed", query)
        print("count:", count)
        _, idlist = search("pubmed", query, count)
        get_abstract("pubmed", idlist)
        merge_txt("download")
    word_dict = build_word_dict("results.txt")
    get_csv(word_dict, "rpb_gene_count.csv")
    wordcloud(word_dict)

if __name__ == '__main__':
    main()


