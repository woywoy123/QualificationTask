from os import walk

likeli = "Normalization"
cutoff = "2600-2800"
Flost = "FLost3"
FLost = {}
FLost["FLost2"] = {}
FLost["FLost3"] = {}

for layer in ["IBL", "Blayer", "layer1", "layer2"]:
    Content = {}
    for (dirpath, dirnames, filenames) in walk("./Output"):
        
        if "FLost" not in dirpath:
            continue
        
        for f in filenames:
            if "Template" not in f or layer not in f or "TRUTH" in f:
                continue
            File = open(dirpath + "/" +f, "r")
            
            Content[dirpath + "/" +f] = []
            for l in File.readlines():
                Content[dirpath + "/" + f].append(l)
    
    
    Truth = {}
    Likelihood = {}
    
    energy = []
    for i in Content:
        Entry = Content[i]
        
        Likelihood[i] = []
        for k in Entry:
            
            if "Truth" in k:
                Truth[i] = k
                continue
            elif "0-" in k and "_GeV" not in k:
                energy = k.split("|")[1:]
                continue
    
            elif "Algorithm" in k or "__" in k or "_GeV" in k or ".txt" in k or k == "\n":
                continue
            
            Likelihood[i].append(k)
    
    

    for t in Truth:
        if "TRUTH" in t:
            continue
        comb = t.split("/")[2]
        FL = t.split("/")[4].split("_")[2].split(".")[0]
        Layer = t.split("/")[4].split("_")[1]
        
        try:
            FLost[FL][Layer]
        except KeyError:
            FLost[FL][Layer] = {}
    
        try: 
            FLost[FL][Layer][comb]
        except KeyError:
            FLost[FL][Layer][comb] = {}
    
        tru = Truth[t].split(" | ")[1:]
        algo = Likelihood[t]
        
        for k in algo:
            pred = k.split(" | ")[1:]
            a = k.split(" | ")[:1][0].replace(" ", "")
             
            i = 0
            abs_er = 0
            for it in range(len(energy)-1):
                tr = tru[it]
                pr = pred[it] 
                ener = energy[it]
                if float(tr) == 0:
                    continuea
                if cutoff in ener:
                    break
                print(ener)
                abs_er += (abs(abs(float(tr))-abs(float(pr)))/float(tr))*100
                i+=1
            FLost[FL][Layer][comb][a] = abs_er/i
    
    Algos = []
    Combi = []
    for L in FLost[Flost]:
        for m in FLost[Flost][L]:
            for algo in FLost[Flost][L][m]:
                if algo not in Algos:
                    Algos.append(algo) 
            if m not in Combi:
                if "TRUTH" in m:
                    continue
                Combi.append(m)
    
    buff = "                        "
    title = ""+buff+" | "
    print("---- Layer ---> ", L)
    trigger = False
    Average = {}
    Minimizer = []
    FitTo = []
    for i in Algos:
        if trigger == False:
            for j in Combi:
                title += " "*(len(j) - len(j))
                title += j
                title += " | "
                Average[j] = 0
    
                if "FitTo" in j:
                    FitTo.append(j)
                if "Minimizer" in j:
                    Minimizer.append(j)
            trigger = True
            print(title) 
        title = " "*(len(buff)-len(i))
        title += i
        title += (" | ")
        
        for j in Combi:
            title += " "*(len(j) - len(str(round(FLost[Flost][layer][j][i], 3))))
            title += str(round(FLost[Flost][layer][j][i], 3))
            title += (" | ")
            Average[j] += round(FLost[Flost][layer][j][i], 3)
        print(title)
    
    print("=======")
    title = ""+buff+" | "
    for i in Average:
        title += " "*(len(i) - len(str(round(Average[i]/len(Algos), 3))))
        title += str(round(Average[i]/len(Algos), 3))
        title += (" | ")
    print(title)
    
    print("")
    print("______")

#for i in Algos:
#    print("--> ", i)
#    best = 0
#    for f, m in zip(FitTo, Minimizer):
#        print(round(FLost[Flost][layer][f][i], 3), round(FLost[Flost][layer][m][i], 3), round(FLost[Flost][layer][f][i] - FLost[Flost][layer][m][i], 3), m.replace("_Minimizer", ""))
#        best += FLost[Flost][layer][f][i] - FLost[Flost][layer][m][i]
#
#    if best > 0:
#        print("FitTo", round(best, 3))
#    else:
#        print("Minimizer", round(best, 3))
#    print("\n")
#
#
#entries_m = []
#entries_f = []
#for i in Algos:
#    print("--> ", i)
#    for f, m in zip(FitTo, Minimizer):
#    	if "To_Smooth" in f:
#    	    r = abs(round(FLost[Flost][layer]["_FitTo"][i] - FLost[Flost][layer][f][i], 3))
#    	    entries_f.append(r)
#
#    	if "r_Smooth" in m:
#       	    r = abs(round(FLost[Flost][layer]["_Minimizer"][i] - FLost[Flost][layer][m][i], 3))
#       	    entries_m.append(r)
#s_f = 0
#s_m = 0
#for i, j in zip(entries_m, entries_f):
#    s_m +=i
#    s_f +=j

#print("Error Increase: ", float(s_f / len(entries_f)), float(s_m / len(entries_m)))


for lay in FLost[Flost]:
    best = 100000
    combi = ""
    Min = ""
    for com in FLost[Flost][lay]:
        for mini in FLost[Flost][lay][com]:
            err = FLost[Flost][lay][com][mini]

            if err < best: 
                best = err
                combi = com
                Min = mini
    print("---"+lay+"-->", Min, combi, best)
