import count_freqs
import sys


def trigramProb(y,y1,y2,trigramLookup,bigramLookup):
    if (y2,y1,y) in trigramLookup and (y2,y1) in bigramLookup:
        return trigramLookup[y2,y1,y]/bigramLookup[y2,y1]
    else:
        return 0
    

def emissionprob(word,tag,pairLookup, unigramLookup, wordLookup):
   
    #e(x|y) = count(y->x)/count(y)  
    if (word,tag) not in pairLookup:
        if isRare(word,wordLookup):
            return float(pairLookup[("_RARE_",tag)]) / float(unigramLookup[tag])
        else:
            return 0
    #print word,tag,pairLookup[(word,tag)],unigramLookup[tag]
    return float(pairLookup[(word,tag)]) / float(unigramLookup[tag])

def initialize(countspath):
    f = open(countspath,"r")
    pairLookup = {}
    wordLookup = {}
    unigramLookup = {}
    bigramLookup = {}
    trigramLookup = {}
    for line in f:
        items = line.replace("\n","").split(' ')
        if items[1] == "WORDTAG":
            assert (items[3],items[2]) not in pairLookup, "Already had that pair!"
            if items[3] in wordLookup:
                wordLookup[items[3]] += float(items[0])
            else:
                wordLookup[items[3]] = float(items[0])
            pairLookup[(items[3],items[2])] = float(items[0])
            # print "(" + items[3],items[2]+ ")","=",items[0]
        else:
            #this is a n-gram
            if items[1][0]=="1":
                #unigram
                unigramLookup[items[2]] = float(items[0])
                #print items[2], "assigned", items[0], "in the unigram lookup"
            elif items[1][0] == "2":
                #bigram
                bigramLookup[(items[2],items[3])] = float(items[0])
            else:
                #trigram
                trigramLookup[(items[2],items[3],items[4])] = float(items[0])
    f.close()

    
    return (pairLookup,wordLookup,unigramLookup,bigramLookup,trigramLookup)

def unigramTagger(word,tags):
    #return sorted(tags, emissionprob(word,tag,pairLookup,unigramLookup,wordLookup))[0]
    maxtag = "O"
    maxpb = 0.0
    curpb = 0.0
    for tag in tags:
        curpb = emissionprob(word,tag,pairLookup,unigramLookup, wordLookup)
        if curpb > maxpb:
            maxpb = curpb
            maxtag = tag
    return maxtag

def devpreprocess(devpath, wordLookup):
     with open(devpath, "r") as f:
        f2 = open(devpath + ".mod","w")
        for line in f:
                    if isRare(line.replace('\n','')):
                        f2.write("_RARE_\n")
                    else:
                        f2.write(line)
        f2.close()

def getsentences(path):
    sents = [[]]
    with open(path,"r") as f:
        for line in f:
            items = line.split(' ')
            if items[0].replace('\n','') == "":
                #end of sent, add STOP and *,*
                sents[-1].append("-STOP-")
                sents[-1].insert(0,"*")
                sents[-1].insert(0,"*")
                sents.append([])
            else:
                sents[-1].append(items[0].replace('\n', ''))
    return sents

def trainpreprocess(trainpath, wordLookup):
    with open(trainpath, "r") as f:
        f2 = open(trainpath + ".mod","w")
        for line in f:
                    if line == '\n':
                        continue
                    items = line.split(' ')
                    if isRare(items[0],wordLookup):
                        f2.write("_RARE_ " + items[1])
                    else:
                        f2.write(line)
        f2.close()

def isRare(word, wordLookup):
    return word not in wordLookup or wordLookup[word] < 5 or word == "_RARE_"

def tagFile(geneDevFile, outfile):
    with open(geneDevFile,"r") as f:
        with open(outfile,"w") as out:
            for line in f.readlines():
                word = line.replace("\n","")
                if word == "":
                    out.write("\n")
                    continue
                out.write(word + " " + unigramTagger(word,tags) + "\n")

def getCounts(intputfile, outputfile):
    i = file(inputfile,"r")
    j = file(outputfile,"w")
    counter = count_freqs.Hmm(3)
    counter.train(i)
    counter.write_counts(j)
    j.close()
    i.close()

# (pairLookup,tagLookup,wordLookup,unigramLookup,bigramLookup,trigramLookup)
getCounts("gene.train","counts.txt")
(pairLookup,wordLookup,unigramLookup,bigramLookup,trigramLookup) = initialize("counts.txt")
trainpreprocess("gene.train",wordLookup)
getCounts("gene.train.mod","counts_mod.txt")
(pairLookup,wordLookup,unigramLookup,bigramLookup,trigramLookup) = initialize("counts_mod.txt")
#print pairLookup[("_RARE_","O")],pairLookup[("_RARE_","I-GENE")]
tags = ["I-GENE","O"]
#devpreprocess("gene.dev",wordLookup)
tagFile("gene.test", "gene_test.p1.out")
#print rarect

