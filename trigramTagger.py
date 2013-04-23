import count_freqs
import sys
import eval_gene_tagger as ev

def trigramProb(y2,y1,y):
    global trigramLookup
    global bigramLookup
    
    if (y2,y1,y) in trigramLookup and (y2,y1) in bigramLookup:
        return trigramLookup[y2,y1,y]/bigramLookup[y2,y1]
    else:
        return 0
    
def emissionProb(word,tag):
    global pairLookup
    global unigramLookup
    global wordLookup
    
    #e(x|y) = count(y->x)/count(y)  
    if (word,tag) not in pairLookup:
        if isRare(word):
            return float(pairLookup[(rareFeat(word),tag)]) / float(unigramLookup[tag])
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


def getsentences(path):
    sents = [[]]
    with open(path,"r") as f:
        for line in f:
            items = line.split(' ')
            if items[0].replace('\n','') == "":
                #end of sent, add STOP and *,*
                #sents[-1].append("-STOP-")
                #sents[-1].insert(0,"*")
                #sents[-1].insert(0,"*")
                sents.append([])
            else:
                sents[-1].append(items[0].replace('\n', ''))
    return sents

def trainpreprocess(trainpath):
    global wordLookup
    
    with open(trainpath, "r") as f:
        f2 = open(trainpath + ".mod","w")
        for line in f:
                    if line == '\n':
                        f2.write('\n')
                        continue
                    items = line.split(' ')
                    if isRare(items[0]):
                        f2.write(rareFeat(items[0]) + ' ' + items[1])
                    else:
                        f2.write(line)
        f2.close()

def isRare(word):
    global wordLookup
    return word not in wordLookup or wordLookup[word] < 5 or word in ['_RARE_','_NUM_','_CAP_','_LASTCAP_']

def rareFeat(word):
    #return the rare feature for this word.
    for c in word:
        if c.isdigit():
            return "_NUM_"
    if word.isupper():
        return "_CAP_"

    if word[-1].isupper():
        return "_LASTCAP_"

    return "_RARE_"

def getCounts(inputfile, outputfile):
    i = file(inputfile,"r")
    j = file(outputfile,"w")
    counter = count_freqs.Hmm(3)
    counter.train(i)
    counter.write_counts(j)
    j.close()
    i.close()

def states(index):
    if index == -1 or index == 0:
        return ['*']
    else:
        return ['I-GENE','O']
    
def tagViterbi(sentence):
    n=len(sentence)
    pi= {(0,'*','*'):1}
    bp= {(0,'*','*'):"O"}
    for k in range(1,n+1):
        s1 = states(k-1)
        s2 = states(k-2)
        s0 = states(k)
        
        for y in s0:
            for y1 in s1:
                #calculate the max and argmax of w in S_k-2
                maxprob, curprob, state = -1,-1,''
                for y2 in s2:
                    curprob = (pi[(k-1,y2,y1)] * trigramProb(y2,y1,y) * emissionProb(sentence[k-1],y))
                    #print "pi[" + str(k-1) + "," + y2 + "," + y1+ "] =",pi[(k-1,y2,y1)]
                    #print "trigramProb(" + y2 + ',' + y2 + ',' + y + ")=", trigramProb(y2,y1,y)
                    #print "e(" + sentence[k-1] + "|" + y + ")=",emissionProb(sentence[k-1],y)
                    #print "---------------------------------------------------"
                    #print "Trigram:", y2,y1,y, "Word:",sentence[k-1]
                    if curprob >= maxprob:
                        state = y2
                        maxprob = curprob
                #print "\n+++\nProb:",maxprob,'\n+++\n'
                
                pi[(k,y1,y)] = maxprob
                bp[(k,y1,y)] = state
                 
    
    
    Sn = states(n)
    Sn1 = states(n-1)
    prevMax = -1
    Tn,Tn1 = '',''
    
    #for k in pi.iterkeys():
    #    print "pi[",k,"] = ",pi[k]
        
    for Yn in Sn:
        for Yn1 in Sn1:
            current = pi[(n,Yn1,Yn)] * trigramProb(Yn1,Yn,'STOP')
            if current >= prevMax:
                Tn = Yn
                Tn1 = Yn1
                prevMax = current
    tags = []
    tags.insert(0,Tn)
    tags.insert(0,Tn1)
    
    for k in range(n-2,0,-1):
        tags.insert(0,bp[(k+2,tags[0],tags[1])])

    #for i in sorted(pi.keys()):
    #    print 'pi',i,pi[i],"argmax =",bp[i]
    #print pi[(3,"I-GENE","I-GENE")]
    return tags

def evalTags(keyFile,predictionsFile):
    gs_iterator = ev.corpus_iterator(file(keyFile))
    pred_iterator = ev.corpus_iterator(file(predictionsFile), with_logprob = False)
    evaluator = ev.Evaluator()
    evaluator.compare(gs_iterator, pred_iterator)
    evaluator.print_scores()

def tagFile(infile, outfile):
    #read the input
    sents = getsentences(infile)
    with open(outfile,"w") as out:
        for sentence in sents:
            #tag and then output.
            tags = tagViterbi(sentence)
            for i in range(0,len(sentence)):
                out.write(sentence[i] + " " + tags[i] + "\n")
            out.write('\n')

    
#(pairLookup,tagLookup,wordLookup,unigramLookup,bigramLookup,trigramLookup)
#getCounts("gene.train","counts.txt")
#(pairLookup,wordLookup,unigramLookup,bigramLookup,trigramLookup) = initialize("counts.txt")
#trainpreprocess("gene.train")
#getCounts("gene.train.mod","counts_mod.txt")
(pairLookup,wordLookup,unigramLookup,bigramLookup,trigramLookup) = initialize("counts_mod.txt")
tagFile("gene.test","gene_test.p3.out")
#tagFile("gene.debug","gene_debug.p2.out")
#evalTags("gene.key","gene_dev.p3.out")
#print pairLookup[("_RARE_","O")],pairLookup[("_RARE_","I-GENE")]
#w = ['STAT5A','mutations','in','the','Src','homology','2','(','SH2',')','and','SH3','domains','did','not','alter','the','BTK','-','mediated','tyrosine','phosphorylation','.','']
#s = tagViterbi(w)
#w.append("End.")
#for i in range(0,len(w)):
#    print w[i],'-',s[i]
#tags = ["I-GENE","O"]
#tagFile("gene.test", "gene_test.p1.out")
#print rarect

