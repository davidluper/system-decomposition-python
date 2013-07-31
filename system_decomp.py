#############################################################################################
#System Decomposition
#by davidluper
#
#Copyright 2010 davidluper
#This file is part of System Decomposition.
#
#System Decomposition is free software: you can redistribute it and/or modify it under the
#terms of the GNU General Public License as published by the Free Software Foundation,
#either version 3 of the License, or (at your option) any later version.
#
#System Decomposition is distributed in the hope that it will be useful, but WITHOUT ANY
#WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#PURPOSE. See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with
#System Decomposition.  If not, see http://www.gnu.org/licenses/.
#############################################################################################

import sys;
import random;
import os;
import shutil;
import SocketServer;
from numpy import *;
from time import strftime;

#TODO
#  DVM SEARCH
#    -  only search chain combinations that have the same input/output nodes
#    -  when making left and right sides of equals sign only enumerate through half the possible combinations (right
#         now every possible combination is put on the left side which is searching through to many, only a specific half need
#         to be put on the left side because the other half appear on the right side)


#CONSTANTS
ENVIRONMENT = "-1";
INGNORE = "#";

MULTI = 0;
ACTUALAVG = 1;
RANDMULTI = 2;
#

directory = "k:\\out\\";
matrixfile = "matrix.txt";
fluxmatrixoutfile = "fluxmatrix.txt"
matrixfile = "matrix.txt";
completedpathsfile = ["paths.dat","paths2.dat"];
grammarfile = "grammar.txt";
fluxlabelfile = "fluxlabel.txt";
fluxfile = "flux.txt";
stoichiometric = int(0);

def TranslateFlux(flux, shift, delim):
    rtn = "";
    for val in flux.split(delim):        
        rtn += str(int(val) + shift) + delim;
        
    return rtn;

def TranslateToAugment(s, delim, g2augmentmap):
    rtn = "";
    trailer = None;    
    for val in s.strip(" ").split(delim):
        if trailer is None or str(g2augmentmap[int(val)]) != trailer:
            rtn += str(g2augmentmap[int(val)]) + delim;
            trailer = str(g2augmentmap[int(val)]);            

    return rtn.strip(delim);

def TranslateFromSymbol(s, delim, symboltointmap):
    rtn = "";
    for val in s.split(delim):        
        rtn += str(symboltointmap[val]) + delim;
    
    return rtn.strip(delim);
    
def LoadStoichiometricMatrix():
    inmatrix = [];

    file = open(directory + matrixfile, 'r');
    for line in file:
        inmatrix.append(line.strip("\r").strip("\n").split(' '));

    #initialize the graph (graph is the adjacency matrix where the first row and first column are the environment)
    graph = [];
    #len + 1 because we need to make room for the environment in the adjacency matrix
    for i in range(len(inmatrix) + 1):
        graph.append([]);

    for i in range(len(inmatrix) + 1):
        for j in range(len(inmatrix) + 1):
            graph[i].append(int(-1));

    #go through the columns of the stoichiometric matrix to translate it into the adjacency matrix
    for col in range(len(inmatrix[0])):
        colsum = 0;
        rowmarker = 0;
        columnmarker = 0;
      
        for row in range(len(inmatrix)):
            val = int(inmatrix[row][col]);
            colsum += val;
            #if val is -1 then it marks a row position in the adjacency matrix
            if (val == -1):
                rowmarker = row + 1;
            #if val is 1 then it marks the column position in the adjacency matrix
            if (val == 1):
                columnmarker = row + 1;

        #mark the edge in the graph
        graph[rowmarker][columnmarker] = col;

    return graph;

def LoadFluxIndexer():
    flux_indexer = {};
    file = open(directory + fluxfile, 'r');
    for line in file:
        line = line.strip("\r").strip("\n");
        flux_indexer[line] = bool(1);

    return flux_indexer;
        
def SaveGrammar(GRAMMAR, X):
    #output the grammar
    cnt = int(0);
    file = open(directory + grammarfile, 'w');
    for item in GRAMMAR:
        s = ""
        for tup in item:
            s+=str(tup[0]) + " " + str(tup[1]) + "|";
        file.write(str(cnt)+">"+s.strip('|')+"\n");
        cnt+=1;
    file.close();
    #for

    file = open(directory + fluxlabelfile, 'w');
    for k,v in X.items():
        file.write(str(k[0]) + "," + str(k[1]) + ":" + str(v[0]) + "," + str(v[1]) + "\n");

    file.close();

def LoadGrammar():
    GRAMMAR = [];
    file = open(directory + grammarfile, 'r');
    for line in file:
        line = line.strip("\r").strip("\n");
        GRAMMAR.append([]);
        rule = line.split('>');
        if len(rule[1].strip(' ')) > 0:
            for s in rule[1].split('|'):
                vals = s.split(' ');
                if vals[1].find("-") > -1:
                    val = vals[1];
                else:
                    val = int(vals[1]);
                    
                GRAMMAR[len(GRAMMAR)-1].append((int(vals[0]),val));

    X = {};
    file = open(directory + fluxlabelfile, 'r');
    for line in file:
        line = line.strip("\r").strip("\n");
        vals1 = line.split(':');
        tup1 = vals1[0].split(',');
        tup2 = vals1[1].split(',');
        X[(int(tup1[0]),int(tup1[1]))] = (int(tup2[0]),int(tup2[1]));
        
    return GRAMMAR, X;

def GetFluxIDFromIndexer(flux, fluxindexer):
    nodelist = flux.split(',');
    length = len(nodelist);
    for i in range(length):
        s = "";
        for j in range(length):                    
            s += "," + nodelist[(i + j) % length];

        if s.strip(",") in fluxindexer: return fluxindexer[s.strip(",")];    

def BuildGrammar(graph, flux_indexer, environ_ndx):
    #build the grammar to use in the CYK algrothim
    counter = int(len(graph) + 1);
    
    #clear grammar
    GRAMMAR = [];
    X = {};
    
    for i in range(counter):
        GRAMMAR.append([]);

    #fluxidcounter = int(-1);
    rtnflux_indexer = {};
    rtnreverseflux_indexer = {};
    for s in flux_indexer:
        #fluxidcounter += 1;
        rtnflux_indexer[s] = flux_indexer[s];#fluxidcounter;
        rtnreverseflux_indexer[flux_indexer[s]] = s; #[fluxidcounter] = s;
        
        flux = TranslateFlux(s, 1, ",").strip(',');        
        rules = [];
        ischain = bool(flux[0] == str(environ_ndx));
        
        if ischain: rules.append(flux.replace(",", " ").strip() + " " + flux[0]);        
        else:
            nodelist = flux.split(',');
            length = len(nodelist);            
            for i in range(length):
                grammaritem = "";
                for j in range(length):                    
                    grammaritem += " " + nodelist[(i + j) % length];

                rules.append(str(grammaritem + " " + grammaritem.strip().split(" ")[0]).strip());
        #else

        for rule in rules:
            nodes = rule.strip().split(" ");
            GRAMMAR[int(nodes[0])].append((int(nodes[0]), counter));
            graphnode = int(nodes[0]);
            start = counter;
            for ndx in range(1, len(nodes) - 2):
                counter += 1;
                GRAMMAR.append([(int(nodes[ndx]), counter)]);                
            #for
            
            X[(graphnode,start)] = (counter,flux_indexer[s]);#fluxidcounter);  
            GRAMMAR.append([(int(nodes[len(nodes) - 2]), int(nodes[len(nodes) - 1]))]);
            counter += 1;            
        #for
    #for

    #make new start state and make a production that turn the new start state into the environment and the enviroment into the start state    
    GRAMMAR[0].append((0,environ_ndx));
    GRAMMAR[environ_ndx].append((environ_ndx,0));
    X[(0,environ_ndx)] = (-1,-1);
    X[(environ_ndx,0)] = (-1,-1);

    #put the literals at the end of the GRAMMAR
    #literals are any node and a - (ex. 0-,1-)
    for node in range(0, len(graph) + 1): GRAMMAR.append([(node, str(node) + "-")]);
    
    return GRAMMAR, X, rtnflux_indexer, rtnreverseflux_indexer;

def ComputeTable(instance, GRAMMAR, literalsndx):
    w = instance.strip().split(" ");
    
    n = len(w);
    r = len(GRAMMAR);
    P = [[[[]]]];
    
    #P = [[[[] for col in range(n)] for row in range(n+1)] for depth in range(r)];
    #init the 3d array of lists
    for col in range(n):
        P.append([[[]]]);
        for row in range(n + 1):
            P[col].append([[]]);
            for depth in range(r):
                P[col][row].append([]);
            #for
        #for
    #for
    
    for i in range(n):
        for j in range(literalsndx, r):
            if GRAMMAR[j][0][1] == w[i]:
                P[i][1][GRAMMAR[j][0][0]].append( ((-1,-1,-1,-1),(-1,-1,-1,-1)) );         
        #for        
    #for
                
    for i in range(2, n+1):   #length of span 1 through len(w)+1
        for j in range(n-i+1):  #start of span, the start index in the w array (0 through len(w))
            for k in range(1,i):    #partition of span must partition into 2 parts
                
                for l in range(literalsndx):  #foreach production the does not have a literal
                    for m in range(len(GRAMMAR[l])):
                        
                        production = GRAMMAR[l][m];
                        if len(P[j][k][production[0]]) > 0 and len(P[j+k][i-k][production[1]]) > 0:
                            #print(((j,i,l),(j,k,production[0]),(j+k,i-k,production[1])));
                            P[j][i][l].append( ((j,k,production[0]),(j+k,i-k,production[1])) );    
                    #for
                #for
            #for
        #for
    #for
                            
    return P;
    
def Interpret(instance, P, graphnodescount, X, multiples):    
    
    w = instance.strip().split(" ");    
    n = len(w);

    rtn = {};
    #if len(P[0][n][1]) > 0: #if the instance was able to be interpreted then will have at least one item. start trace in that cell
    if len(P[0][n][0]) > 0: #if the instance was able to be interpreted then will have at least one item. start trace in that cell    
        branchmap = {};        
        queue = [];
        #queue.append( ( [((0,n,1),0,-1)], [], {}) );
        queue.append( ( [((0,n,0),0,-1)], [], {}) );
        while len(queue) > 0:
            stack = queue[0][0];
            interp = queue[0][1]
            visited = queue[0][2];
           
            while len(stack) > 0:                
                node = stack[len(stack) - 1][0];
                height = stack[len(stack) - 1][1];
                
                if P[node[0]][node[1]][node[2]][0][0][0] == -1: #is leaf
                    pos = len(stack)-1;
                    start = stack[pos][0][2];
                    hold = int(-1);
                    #lastheight = stack[pos][1];            
                    graphnode = start;
                    goneleft = bool(0);
                    #while not at root and either hold == -1 (meaning we have not run an initial pass) or are not at a literal or we have not gone left yet
                    #once we are at a literal and have gone left we break. 
                    while pos > 0 and (hold == -1 or graphnode > graphnodescount or not goneleft): #(graphnode == start and not goneleft)):
                        lastpos = pos;  #lastpos is child position on stack
                        while stack[pos][1] == stack[lastpos][1]:
                            pos -= 1;                        
                        #pos is now the parent position on the stack

                        #if the child was a left child then we mark that we have gone left, if we have gone left 
                        if stack[lastpos][2] == 0: goneleft = bool(1);
                            
                        hold = graphnode;
                        graphnode = stack[pos][0][2];

                    #if where we stopped and where we started are the same literal then hold equals graphnodes left child
                    #this will give us the proper flux label
                    if start == graphnode:
                        #print(("*", graphnode));
                        hold = stack[pos+1][0][2];

                    #print (str(X[(graphnode,hold)][1]));
                    interp.append((node[2],str(X[(graphnode,hold)][1])));                    
                    stack.pop();
                else:
                    #make a copy of the current state of the stack for any branch of the reconstruction
                    #then add the branch to the new stack for evaluation later                                   
                    if multiples and (node[0],node[1],node[2],height) not in branchmap:
                        branchmap[(node[0],node[1],node[2],height)] = bool(1);
                        for i in range(1,len(P[node[0]][node[1]][node[2]])):                                 
                            queue.append(([],[],{}));
                            for j in stack:
                                queue[len(queue)-1][0].append(j);
                            for j in interp:
                                queue[len(queue)-1][1].append(j);
                            for _k,_v in visited.items():
                                queue[len(queue)-1][2][_k] = _v;
                                    
                            if (node[0],node[1],node[2],height) not in queue[len(queue)-1][2]: 
                                queue[len(queue)-1][0].append((P[node[0]][node[1]][node[2]][i][1],height+1,1)); 
                                queue[len(queue)-1][0].append((P[node[0]][node[1]][node[2]][i][0],height+1,0));                                
                                queue[len(queue)-1][2][(node[0],node[1],node[2],height)] = bool(1);
                                
                    #current branch will always be at P[n][n][r][0]
                    if (node[0],node[1],node[2],height) not in visited:
                        stack.append((P[node[0]][node[1]][node[2]][0][1],height+1,1));                           
                        stack.append((P[node[0]][node[1]][node[2]][0][0],height+1,0));
                        visited[(node[0],node[1],node[2],height)] = bool(1);
                    else:
                        stack.pop();
            
            if len(interp) == n:          
                s = "";
                for item in interp:
                    s += str(item[1]) + " ";

                s = s.strip(' ');

                if s not in rtn:
                    rtn[s] = bool(1);
                    
            queue.pop(0);

    rtnlist = [];
    for item in rtn:
        rtnlist.append(item);
        
    return rtnlist;
#Interpret
                
#makes a string from a list of items
def MakeString(items, start, count):
    rtn = "";
    for i in range(int(start), int(start) + int(count)):
        if int(i) > int(start):
            rtn += "," + str(items[int(i)]);
        #if
        else:
            rtn += str(items[int(i)]);
        #else
    #for
    return rtn;
#MakeString

def Decompose(graph, environ_ndx, symbolmap, strmap, intmap):
    if strmap == None: strmap = {};
    if intmap == None: intmap = {};
    
    if graph is None:
        if stoichiometric == 1:
            graph = LoadStiociometricMatrix();
        #end stoichiometric translation
        elif stoichiometric == 0:
            #the adjacency matrix needs to already have the environment as row and column 0
            graph,strmap,intmap = LoadAdjacencyMatrix(symbolmap);
    
    queue = [];
    flux_indexer = { };

    #put every compartment in a path of length one and then on the queue
    for node in range(len(graph)):
        starter = [ int(node) ];
        queue.append(starter);
    #for

    #while queue is not empty
    while len(queue) > 0:
        #get the next list of nodes to evaluate from the queue
        node_list = queue.pop(0);
        
        row = int(node_list[len(node_list) - 1]);
        column = int(node_list[0]);
        #true if the last node has an edge to the first node
        
        found = graph[row][column] > -1;
            
        if found:
            #see if this node list has already been found
            exists = bool(0);
            #look for every valid combination of the node list
            #i.e. 1,2,3 - 3,1,2 - 2,3,1
            for i in range(len(node_list)):
                templist = [];
                #make the next combination
                for j in range(i, len(node_list) + i):
                    templist.append(node_list[j % len(node_list)]);
                #for

                #turn the combination into a string
                s = MakeString(templist, 0, len(templist));

                #if this representation of the node list is in the indexer
                #then we have already found this flux
                if s in flux_indexer:
                    exists = bool(1);
                    break;
                #if
            #for

            #if the flux does not exist yet, put it in the indexer
            if not exists:                
                insert = MakeString(node_list, 0, len(node_list));                
                flux_indexer[insert] = int(-1);
            #if
        #if

        #get the edge list for the last node in the node list
        edges = graph[row];

        #iterate through all the edges to look for connected nodes
        for i in range(len(edges)):
            #if there is an edge and the connecting node is not in the node list already
            if edges[int(i)] > -1 and int(i) not in set(node_list):
                #copy the list to a new list
                new_list = [];
                for j in range(len(node_list)):
                    new_list.append(node_list[j]);
                #for           

                #add the next node
                new_list.append(i);
                #put it on the queue
                queue.append(new_list);
            #if
        #for
    #while


    #output fluxes in a text file where each line is a series of nodes representing a flux
    #output the fluxes
    
    #for
    

    #output the flux matrix
    #translate decomposition into flux matrix
    fluxmatrix = [];

    edgecnt = 0;
    edgematrix = [];
    for i in range(len(graph)):
        edgematrix.append([]);
        for j in range(len(graph)):
            if graph[i][j] == 1:
                edgematrix[i].append(edgecnt);
                edgecnt += 1;
            else: edgematrix[i].append(-1);                    
    
    #initialize flux matrix
    for i in range(edgecnt):
        fluxmatrix.append([]);
        for j in range(len(flux_indexer)):
            fluxmatrix[i].append(0);

    #loop through each decomposed flux.  lookup each edge label in the graph
    #and for the particular row we are on set that label to "1"
    
    col = int(0);
    out_chains = [];
    out_cycles = [];
    set_chains = [];
    set_cycles = [];
    for flux in flux_indexer:            
        compartments = flux.split(',')
        edgeset = [];
        for i in range(len(compartments)):
            edgeset.append(edgematrix[int(compartments[i])][int(compartments[(i + 1) % len(compartments)])]);
            row = edgematrix[int(compartments[i])][int(compartments[(i + 1) % len(compartments)])];
            fluxmatrix[row][col] = 1;
            edgeset.sort();
        if int(compartments[0]) == environ_ndx:
            set_chains.append(edgeset);
            out_chains.append(flux);
        else:
            set_cycles.append(edgeset);
            out_cycles.append(flux);
            
        col += 1;

    file = open(directory + fluxfile, 'w');
    for x in range(len(out_chains)):
        flux_indexer[out_chains[x]] = x; #this is important, it aligns the flux number ids in the histogram with the dvm
        file.write(out_chains[x] + "\n");
    for x in range(len(out_cycles)):
        flux_indexer[out_cycles[x]] = len(out_chains) + x; #this is important, it aligns the flux number ids in the histogram with the dvm
        file.write(out_cycles[x] + "\n");
    file.close();
    
    #output the flux matrix
    file = open(directory + fluxmatrixoutfile, 'w');
    for flux in fluxmatrix:
        for item in flux: file.write(str(item) + " ");
        file.write("\n");
    #for
    file.close();

    return graph, flux_indexer, fluxmatrix, set_chains, set_cycles, out_chains, out_cycles, strmap, intmap;

def ListHash(_list):
    rtn = "";
    for x in range(len(_list)):
        if _list[x] != 0.0:
            rtn += str(_list[x]) + "," + str(x) + ":";
    return rtn;

def ComputeCoefficientsAndDVM(flux_indexer,GRAMMAR,X,graph,histogram,startcnt,dvm,_pathfile,_symboltointmap,_inttosymbolmap,g2augmentmap,augment2gmap,addenviron,pathdict,numberofinstances):

    print ("START",strftime("%Y-%m-%d %H:%M:%S"));
    uniquecnt = 0;
    vectorlen = len(flux_indexer);
    #histogram = [0 for i in range(vectorlen)];
    rank = 0;
    if len(dvm) > 0: rank = matrixrank(array(dvm));
    
    #--------------------------------------------
    fluxlens = [0 for i in range(vectorlen)];
    for k,v in X.items():        
        if fluxlens[v[1]] != 0: continue;
        if k[0] == 0: continue;
        #1 is the start graphnode
        #(i[1][0] - i[0][1]) + 1 is the number of elements in the rest of the flux
        length = 1 + ((v[0] - k[1]) + 1);
        #dont do this it throws off the frequency count
        #if chain count the last 1 in the flux as well
        #if k[0] == 1: length += 1;
        fluxlens[v[1]] = length;
    #-------------------------------------------
        
    file = open(directory + _pathfile, 'r');
    cnt = 0;
    dvmdict = {};
    for dv in dvm: dvmdict[MakeString(dv, 0, len(dv))] = bool(1);

    breakcnt = int(0);
    
    for line in file:
        breakcnt += 1;

        if breakcnt > numberofinstances:
            break;
        
        if cnt < startcnt:
            cnt += 1;
            continue;
        
        line = line.strip("\r").strip("\n").strip(' ');                
        interpretation = [];

        if line in pathdict:
            histogram = addvec(histogram,pathdict[line]);
        else:            
            #instance = TranslateFlux("-1 0 1 3 4 1 3 1 2 3 4 1 3 1 2 3 -1", 2, " ").replace(" ", "- ");
            #instance = TranslateFlux("-1 0 1 2 3 4 1 3 1 2 3 -1", 2, " ").replace(" ", "- ");
            #instance = TranslateFlux("-1 0 -1", 2, " ").replace(" ", "- ");            
            #instance = TranslateFlux("-1 2 -1", 2, " ").replace(" ", "- ");

            #we need to translate from symbols into adjacency matrix indexes, in the oyster reef or other PTN data models this usually equates to a shift of + 1
            #then we build the grammar where we add an additional start state so the environment becomes 0 and 1,
            #so 0 in the data set becomes 2 in the grammar, etc and we need to shift + 1 from the adjacency matrix indexes
            instance = (str(TranslateFromSymbol(line, " ", _symboltointmap)), str(TranslateFromSymbol("-1 " + line + " -1", " ", _symboltointmap)))[addenviron];

            #now the symbol ints in instance map to indexes in the origonal graph adjacency matrix
            #if the origonal to augmented adjacency map is not null then translate the instance into something the new grammar will be able to parse
            if g2augmentmap is not None: instance = TranslateToAugment(instance, " ", g2augmentmap);
            
            #perform the final shift to accomadate the extra start symbol contained in the grammar and add the start and stop states at beginning and end            
            instance = str("0 " + TranslateFlux(instance, 1, " ") + "0 ").replace(' ', '- ');
            
            instancelen = len(instance.split(' '));
            if instancelen > _MAXLEN:
                print ("over max",instancelen);
                continue;

            uniquecnt += 1;
        
            #instance = "0- 1- 2- 1- 2- 3- 1- 0-";
            P = ComputeTable(instance, GRAMMAR, len(GRAMMAR) - (len(graph) + 1));
            interpvecs = [];
            interps = [];

            interpdict = {};
            for item in Interpret(instance, P, len(graph),X,bool(histmethod == ACTUALAVG or histmethod == RANDMULTI)):
                tempvec = tointerpvector(item.split(' '), vectorlen, fluxlens);
                s = ListHash(tempvec);
                
                if s not in interpdict:
                    interpdict[s] = bool(1);
                    item = item[3:len(item)-6]; #trim the -1 at the start and the two -1's at the end                    
                    interps.append(item);
                    interpvecs.append(tempvec);

            #any path should be able to be interpreted so this executing means something in the setup was wrong
            if len(interpvecs) < 1: print ("NO INTERPRETATION, " + instance,TranslateFlux("-1 " + line + " -1", 2, " "));
            #if len(interpvecs) > 1: print ("MULTIPLE INTERPRETATIONS");
            
            if histmethod == MULTI:
                histogram = addvec(histogram,interpvecs[0]);
                pathdict[line] = interpvecs[0];
            elif histmethod == ACTUALAVG:
                temphistogram = [0.0 for i in range(len(histogram))];
                for i in range(len(interpvecs)):
                    temphistogram = addvec(temphistogram,scalevec(interpvecs[i], len(interpvecs)));

                pathdict[line] = temphistogram;
                histogram = addvec(histogram,temphistogram);
            elif histmethod == RANDMULTI:
                if len(interpvecs) > 1:
                    z = random.randint(0, len(interpvecs)-1);
                    histogram = addvec(histogram,interpvecs[z]);
                    pathdict[line] = interpvecs[z];                    
                else:
                    histogram = addvec(histogram,interpvecs[0]);
                    pathdict[line] = interpvecs[0];

            if dvm is not None and computeinlinedvm:
                for x in range(len(interpvecs)):
                    for i in range(len(interpvecs)):
                        if x == i:
                            continue;
                        
                        dv = subvec(interpvecs[x],interpvecs[i]);
                        dvcheck = MakeString(dv, 0, len(dv));
                        if dvcheck not in dvmdict:
                            dvmdict[dvcheck] = bool(1);
                            dvm.append(dv);
                            newrank = matrixrank(array(dvm));
                            if newrank > rank: rank = newrank; print "dvm";
                            else: dvm.pop();                 
        #else
                        
        cnt += 1;
        if cnt % 100 == 0: print (cnt,strftime("%Y-%m-%d %H:%M:%S"));

        if cnt % 25 == 0:
            try: os.remove(directory + "log.txt");
            except: t = 0;
            logfile = open(directory + "log.txt", 'w');
            logfile.write(str(cnt));
            logfile.write("\n");

            strhist = ""
            for item in histogram:
                strhist = strhist + (str(item) + ",");
            logfile.write(strhist.strip(',') + "\n");
            
            for item in dvm:
                strdvm = ""
                for x in item:
                    strdvm = strdvm + str(x) + ",";

                logfile.write(strdvm.strip(',') + "\n");

            logfile.close();
    #for
    file.close();

    print ("STOP",strftime("%Y-%m-%d %H:%M:%S"));
    print ("UNIQUE CNT : " + str(uniquecnt));
    return histogram,dvm,pathdict;

def scalevec(x, factor):
    for i in range(len(x)):
        x[i] = float(x[i]) / float(factor);
    return x;

def subvec(x, y):    
    rtn = [0.0 for i in range(len(x))];
    for i in range(len(x)): rtn[i] = float(x[i]) - float(y[i]);
    return rtn;

def addvec(x, y):
    rtn = [0.0 for i in range(len(x))];
    for i in range(len(x)): rtn[i] = float(x[i]) + float(y[i]);
    return rtn;

def tointerpvector(interpretation, vectorlen, fluxlens):
    rtn = [0.0 for i in range(vectorlen)];    
    for i in range(1, len(interpretation) - 2): rtn[int(interpretation[i])] = float(rtn[int(interpretation[i])]+1);    
    for i in range(len(rtn)): rtn[i] = float(rtn[i] / float(fluxlens[i]));
    return rtn;

def matrixrank(A,tol=1e-8):
    s = linalg.svd(A,compute_uv=0)
    return sum( where( s>tol, 1, 0 ) )

def LoadAdjacencyMatrix(symbolmap):    
    intmap = {};
    strmap = {};
    translist = [];
    vals = symbolmap;
        
    if _ADDENVIRON:
        translist.append(ENVIRONMENT);
        intmap[ENVIRONMENT] = int(len(translist)-1);
        strmap[int(len(translist)-1)] = ENVIRONMENT;
        

    for item in vals:
        translist.append(item);
        intmap[item] = int(len(translist)-1);
        strmap[int(len(translist)-1)] = item;

    inmatrix = [];

    file = open(directory + matrixfile, 'r');
    for line in file:
        inmatrix.append(line.strip("\r").strip("\n").strip(" ").split(' '));

    #initialize the graph (graph is the adjacency matrix where the first row and first column are the environment)   
    graph = [];
    for i in range(len(inmatrix)):
        graph.append([]);
            
        for j in range(len(inmatrix[i])):
            if inmatrix[i][j] == "0":
                graph[i].append(-1);
            else:
                graph[i].append(int(inmatrix[i][j]));
                
    return graph, strmap, intmap;

def initAdjacencyMatrix(symbolmap, addenviron):
    transdict = {};    
    translist = [];
    
    intmap = {};
    strmap = {};

    if symbolmap is not None:
        #symbolfile = open(strdir + "symbolmapin.txt", 'r');
        #vals = symbolfile.readline().strip("\r").strip("\n").split(",");
        #symbolfile.close();

        vals = symbolmap;
        
        if addenviron:
            translist.append(ENVIRONMENT);
            intmap[ENVIRONMENT] = int(len(translist)-1);
            strmap[int(len(translist)-1)] = ENVIRONMENT;
            transdict[ENVIRONMENT] = {};

        for item in vals:
            translist.append(item);
            intmap[item] = int(len(translist)-1);
            strmap[int(len(translist)-1)] = item;
            transdict[item] = {};

    return transdict, translist, intmap, strmap;

def buildAdjacencyMatrix(instance, addenviron, symbolmap):
    transdict, translist, intmap, strmap = initAdjacencyMatrix(symbolmap, addenviron);
    _s = instance.strip("\r").strip("\n").strip(" ");
    if addenviron: _s = ENVIRONMENT + " " + _s + " " + ENVIRONMENT;
        
    vals = _s.split(' ');   
       
    for i in range(1, len(vals)):
        if vals[i-1] not in transdict:
            translist.append(vals[i-1]);
            intmap[vals[i-1]] = int(len(translist)-1);
            strmap[int(len(translist)-1)] = vals[i-1];
            transdict[vals[i-1]] = {};

        if vals[i] not in transdict[vals[i-1]]:
            transdict[vals[i-1]][vals[i]] = int(0);
              
        transdict[vals[i-1]][vals[i]] += 1;

    if vals[len(vals)-1] not in transdict:
        translist.append(vals[len(vals)-1]);

    return CompileAdjacencyMatrix(transdict, translist, intmap, strmap, None);

def buildAdjacencyMatrixFromFile(strdir, strfilename, addenviron, symbolmap):
    transdict, translist, intmap, strmap = initAdjacencyMatrix(symbolmap, addenviron);
    
    pathfile = open(strdir + strfilename, 'r');    

    for line in pathfile:
        _s = line.strip("\r").strip("\n").strip(" ");
        if addenviron: _s = ENVIRONMENT + " " + _s + " " + ENVIRONMENT;
        
        vals = _s.split(' ');   
        
        for i in range(1, len(vals)):
            if vals[i-1] not in transdict:
                translist.append(vals[i-1]);
                intmap[vals[i-1]] = int(len(translist)-1);
                strmap[int(len(translist)-1)] = vals[i-1];
                transdict[vals[i-1]] = {};

            if vals[i] not in transdict[vals[i-1]]:
                transdict[vals[i-1]][vals[i]] = int(0);
                
            transdict[vals[i-1]][vals[i]] += 1;

        if vals[len(vals)-1] not in transdict:
            translist.append(vals[len(vals)-1]);
                                  
    pathfile.close();

    return CompileAdjacencyMatrix(transdict, translist, intmap, strmap, strdir);

def CompileAdjacencyMatrix(transdict, translist, intmap, strmap, strdir):
    outmatrix = [];
    outmatrix_weighted = [];
    
    for item in translist:
        outmatrix.append([]);
        outmatrix_weighted.append([]);
        for i in range(len(translist)):
            outmatrix[intmap[item]].append( (int(-1),int(1))[strmap[i] in transdict[item]] );
            if strmap[i] in transdict[item]:
                outmatrix_weighted[intmap[item]].append( transdict[item][strmap[i]] );
            else:
                 outmatrix_weighted[intmap[item]].append( int(0) );            

    if strdir is not None:                                               
        matrixfile = open(strdir + "builtmatrix.txt", 'w');
        matrixfile_weighted = open(strdir + "builtmatrix_weighted.txt", 'w');
        symbolmap = open(strdir + "symbolmapout.txt", 'w');
        cnt = int(0);

        for _list in outmatrix_weighted:
            for _item in _list:
                matrixfile_weighted.write(str(_item) + " ");
            matrixfile_weighted.write("\n");
            
        for _list in outmatrix:
            for _item in _list:
                matrixfile.write(str(_item) + " ");
            matrixfile.write("\n");
            
            symbolmap.write(str(cnt) + "=" + strmap[cnt]);
            symbolmap.write("\n");
            cnt += 1;

        matrixfile_weighted.close();
        matrixfile.close();
        symbolmap.close();

    return outmatrix, outmatrix_weighted, strmap, intmap;

def augmentAdjacencyMatrix(matrix, inlist):
    outmatrix = [[-1 for i in range(len(inlist) + 1)] for j in range(len(inlist) + 1)];
    oldnew = {};
    newold = {};
    
    cnt = int(1);
    for i in range(len(matrix)):
        if i in inlist:
            oldnew[int(i)] = cnt;
            newold[cnt] = int(i);
            cnt += 1;
        else:
            oldnew[int(i)] = 0;        
    
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):            
            row = (0, oldnew[i])[i in inlist];
            col = (0, oldnew[j])[j in inlist];
            if i != j: outmatrix[row][col] = (matrix[i][j], 1)[outmatrix[row][col] == 1];
            else: outmatrix[row][col] = -1;   
    
    return outmatrix, oldnew, newold;

def ReverseMap(augment2gmap, s, delim):
    rtn = "";
    for item in s.split(delim):
        rtn += (delim + str(augment2gmap[int(s)]));
    return rtn;

def PrintMatrix(matrix):
    for i in range(len(matrix)):
        s = "";
        for j in range(len(matrix[i])):
            s += str(matrix[i][j]) + " ";

        print (s + "\n");

    return;

def MakeDVMRow(_len, set1, set2):
    rtn = [int(0) for i in range(_len)];
    for i in range(len(set1)): rtn[set1[i]] = 1;
    for i in range(len(set2)): rtn[set2[i]] = -1;
    
    return rtn;

def MergeEdgeSets(setlist, indexlist):
    rtnlist = [];
    try:
        
        #print (len(setlist));
        #print (indexlist);
        for ndx in indexlist:
            for item in setlist[ndx]:
                rtnlist.append(item);
        
        return rtnlist;
    except:
        print(("ERROR: ", rtndict, len(setlist),indexlist));

def DVMSet(set1, set2):
    if len(set1) != len(set2): return bool(0);

    set1.sort();
    set2.sort();
    for i in range(len(set1)):
        if set1[i] != set2[i]: return bool(0);
    return bool(1);

def Next(_list, mins, maxs):    
    rtnbool = bool(1);
    for i in range(len(_list)):
        pos = len(_list) - (i+1);
        _list[pos] = (_list[pos] + 1) % maxs[pos];
        if _list[pos] != 0:            
            for j in range(pos+1, len(_list)):
                _list[j] = max(_list[pos] + (j-(pos+1)) + 1, mins[j]);
                if _list[j] >= maxs[j]: rtnbool = bool(0);
            break;
        else: _list[pos] = mins[pos];
          
        if pos == 0: rtnbool = (bool(0),bool(1))[_list[0] != mins[0]];
            
    return rtnbool, _list, mins;

def f(l,n):
    rtn = [];
    for item in l: rtn.append(item + n);
    return rtn;

def SplitCycles(lhc, cycles):
    left = [];
    right = [];

    mark = 0;
    for i in range(len(cycles)):
        if mark < len(lhc) and i == lhc[mark]:
            left.append(cycles[i]);
            mark += 1;
        else:
            right.append(cycles[i]);
            
    return left, right;

def DVMSearchReal(set_chains, set_cycles, adj_len, flux_chains, flux_cycles):
    print ("START",strftime("%Y-%m-%d %H:%M:%S"));
    dvm = [];
    rank = -1;
    masterset = [];
    for item in set_chains: masterset.append(item);
    for item in set_cycles: masterset.append(item);

    #print (masterset);
    cycleoffset = len(set_chains);       

    cyclecnt_max = len(set_cycles);

    exitcnt = 0;
    cnt = 0;
    for x in range(2, cyclecnt_max):

        chains = None;  
        chainsmaxs = None;
        chainmins = None;
    
        partition = 1;    
        cnt = 0;
        chaincnt = 0;
        done = bool(0);

        print(("CHECK " + str(x), strftime("%Y-%m-%d %H:%M:%S")));
        while not done:
            cycles = [];# [0, 1];    
            cyclemaxs = [];# [len(set_cycles)-1, len(set_cycles)];
            cyclemins = [];
            for i in range(0, x):
                cycles.append(i);
                cyclemins.append(i);
                cyclemaxs.append(len(set_cycles)-(x-(i+1)));

            lhc = []; 
            lhcmaxs = [];
            lhcmins = [];

            for i in range(0, partition):
                lhc.append(i);
                lhcmins.append(i);
                lhcmaxs.append(len(cycles)-(partition-(i+1)));
                
            lhcmaxs[0] = lhcmaxs[len(lhcmaxs)-1] / 2;
            
            while 1:#enumerate through each combination of unique cycle set of length 2 to n (n is total number of cycles)
                checklist = f(cycles, cycleoffset);
                if chains is not None:
                    checklist += chains;

                checksetlist = MergeEdgeSets(masterset, checklist);
                checkset = bool(1);

                checkdict = {};
                for item in checksetlist:
                    if item in checkdict:
                        checkdict[item] += 1;
                    else:
                        checkdict[item] = int(1);

                for item in checkdict:
                    if checkdict[item] % 2 != 0:
                        checkset = bool(0);
                        break;
                    
                #print (checkdict);
                #sys.exit();

                if checkset:
                    while 1:#make left and right sides of equals sign
                        left = None;
                        right = None;

                        left, right = SplitCycles(lhc, cycles);
                        left = f(left, cycleoffset);
                        right = f(right, cycleoffset);
                        
                        if chains is not None:
                            left = chains[0:1] + left;
                            right = chains[1:2] + right;

                        #analyze and compare
                        if DVMSet(MergeEdgeSets(masterset, left), MergeEdgeSets(masterset, right)):                    
                            dv = MakeDVMRow(len(masterset), left, right);

                            dvm.append(dv);
                            newrank = matrixrank(array(dvm));
                            if newrank > rank:
                                #print (dv);
                                rank = newrank;
                                print (("dv", left, right));
                                linel = "";
                                liner = "";
                                for item in left: linel += str(item) + ",";
                                for item in right: liner += str(item) + ",";
                                _file = open(directory + "dv-log.txt", 'a');
                                _file.write(linel + " : " + liner);
                                _file.write("\n");
                                _file.close();
                            else: dvm.pop();

                        rtn1, lhc, lhcmins = Next(lhc, lhcmins, lhcmaxs);
                        if not rtn1:                        
                            break;
                    
                #get next
                rtn, cycles, cyclemins = Next(cycles, cyclemins, cyclemaxs);
                                    
                #print((cycles,cyclemins, cyclemaxs));
                if not rtn:
                    #print(("CYCLE DONE",x,partition,strftime("%Y-%m-%d %H:%M:%S")));
                    cnt += 1;                   
                    partition += 1;
                    if partition > x / 2:
                        partition = 1;

                        lhc = [];
                        lhcmaxs = [];
                        lhcmins = [];

                        for i in range(0, partition):
                            lhc.append(i);
                            lhcmins.append(i);
                            lhcmaxs.append(len(set_cycles)-(partition-(i+1)));
                
                        if chains == None:
                            chains = [0, 1];    
                            chainsmaxs = [len(set_chains)-1, len(set_chains)];
                            chainmins = [0, 1];
                            rtn = bool(1);
                        else:
                            rtn, chains, chainmins = Next(chains, chainmins, chainsmaxs);

                            while 1:
                                _0start = str(flux_chains[chains[0]][0:1]);
                                _1start = str(flux_chains[chains[1]][0:1]);
                                _0end = str(flux_chains[chains[0]][len(flux_chains[chains[0]]) - 1:len(flux_chains[chains[0]])]);
                                _1end = str(flux_chains[chains[1]][len(flux_chains[chains[1]]) - 1:len(flux_chains[chains[1]])]);
                                if _0start != _1start or _0end != _1end:
                                    rtn, chains, chainmins = Next(chains, chainmins, chainsmaxs);
                                    if not rtn:
                                        break;
                                else:
                                    break;
                                    
                            _file2 = open(directory + "chainlog.txt", 'a');
                            _file2.write(str(chains[0]) + " : " + str(chains[1]));
                            _file2.write("\n");
                            _file2.close();
                                
                        chaincnt += 1;
                            
                        #print(("NEXT CHAIN",chaincnt,strftime("%Y-%m-%d %H:%M:%S")));

                        if not rtn: done = bool(1);                    
                        
                    break;                           
    
    print ("STOP",strftime("%Y-%m-%d %H:%M:%S"));
    print("number of possible cycle combinations: " + str(cnt));
    print("difference vectors found: " + str(len(dvm)));
    
    return dvm;

def MakeDataSet(path, datafile, numinstances):
    pathways = [];

    _file = open(datafile, 'r');
    
    for line in _file:
        pathways.append(line);
        
    _file.close();

    mark = int(float(len(pathways)) * 0.25);

    _dict = {};

    while len(_dict) < numinstances:
        ndx = random.randint(mark, len(pathways)-1);
        if int(len(pathways[ndx].split(' '))) < _MAXLEN:
            _dict[ndx] = bool(1);
        
    _file = open(path + "paths.dat", 'w');
    _file2 = open(path + "stats.csv", 'w');
    
    for item in _dict.keys():
        _file.write(pathways[int(item)]);
        _file2.write(str(item) + "," + str(len(pathways[int(item)].split(' '))) + "\n");
        
    _file.close();
    _file2.close();
    return;

def LoadHistogram(histfilepath):
    rtn = [];#[0.0 for i in range(vectorlen)];
    histfile = open(histfilepath, 'r');
    
    for line in histfile:
        val = line.strip("\r").strip("\n");
        rtn.append(float(val));
        
    histfile.close();
    #print rtn;
    return rtn;

def DVMDistanceCalc(a, ahat, V):    
    s = matrix(V * V.T).I * V * matrix(a - ahat).T;
    p = ahat + (s.T * V);
    return EuclideanDistance(a, p), p;
    #return a - p, p;
def EuclideanDistance(a, ahat):
    dist = float(0);
    for i in range(a.shape[1]):
        dist += ((a[0, i] - ahat[0, i]) ** 2.0);             
    return dist ** .5;

def LoadInlineDVM (dvmpath):
    dvmfile = open(dvmpath, 'r');
    dvmlist = [];
    for line in dvmfile:
        dvmlist.append(line.strip().split(','));
        
    dvmfile.close();

    rtn = matrix(zeros(shape=(len(dvmlist), len(dvmlist[0]))));
    
    for i in range(len(dvmlist)):
        for j in range(len(dvmlist[i])):
            rtn[i, j] = float(dvmlist[i][j]);

    return rtn;
    
def LoadDVM(dvmpath, fluxpath):
    fluxfile = open(fluxpath);
    cnt = 0;
    for line in fluxfile:
        cnt += 1;

    fluxfile.close();
    
    dvmfile = open(dvmpath, 'r');
    dvmlist = [];
    for line in dvmfile:
        vals = line.strip().split(':');
        if len(vals) is not 2: continue;
        vl = vals[0].strip(" ").strip(",").split(',');
        vr = vals[1].strip(" ").strip(",").split(',');
        dvmlist.append((vl, vr));
        
    dvmfile.close();
    rtn = matrix(zeros(shape=(len(dvmlist), cnt)));
    
    for i in range(len(dvmlist)):
        for ndx in dvmlist[i][0]: rtn[i, int(ndx)] = 1.0;
        for ndx in dvmlist[i][1]: rtn[i, int(ndx)] = -1.0;

    return rtn;

def ReverseTransform(strdir, adjmatrix, keepfluxes):
    freqfile = open(strdir + "freqcount.txt", 'r');
    freqfluxfile = open(strdir + "freqcount_flux.txt", 'r');

    #for line in freqfile:
    cnt = int(0);
    while 1:
        val = freqfile.readline().strip("\n");
        
        if val == "":
            return adjmatrix;

        edges = freqfluxfile.readline().strip("\n").split(',');

        if cnt not in keepfluxes:
            cnt += 1;
            continue;
        
        for i in range(len(edges)):
            row = int(edges[i]);
            col = int(edges[(i+1)%len(edges)]);
            adjmatrix[row][col] += float(val);

        cnt += 1;
    return adjmatrix;

def LoadAndZeroAdj(strdir, matrixfilename):
    matrixfile = open(strdir + matrixfilename, 'r');
    rtn = [];
    cnt = 0;
    for line in matrixfile:
        rtn.append([]);
        
        for item in line.split(' '):            
            if item.strip() != "": rtn[cnt].append(float(0.0));
        
        cnt += 1;

    matrixfile.close();

    return rtn;
#----------------------------------------------------------------------------------------------------
         
random.seed(None);

decompose = bool(0);
loadgraph = bool(0);
loadfluxindexer = bool(0);
buildgrammar = bool(0);
loadgrammar = bool(0);
savegrammar = bool(0);
coefficients = bool(0);
histmethod = int(0);
_MAXLEN = 100;
dmine = bool(0);
_ADDENVIRON = bool(0);
loadlog = bool(0);
symbol_map_= None;
DVM = bool(0);
instancenum = int(999999);
makeadj = bool(0);
computeinlinedvm = bool(0);

for x in range(len(sys.argv)):
    if sys.argv[x] == "-dvm":
        DVM = bool(1);
        decompose = bool(1);
        makeadj = bool(1);
        _ADDENVIRON = bool(1);
    if sys.argv[x] == "-decompavg":
        decompose = bool(1);
        buildgrammar = bool(1);
        coefficients = bool(1);
        makeadj = bool(1);
        histmethod = ACTUALAVG;
        _ADDENVIRON = bool(1);
    if sys.argv[x] == "-decompinlinedvm":
        decompose = bool(1);
        buildgrammar = bool(1);
        coefficients = bool(1);
        makeadj = bool(1);
        histmethod = ACTUALAVG;
        _ADDENVIRON = bool(1);
        computeinlinedvm = bool(1);
    if sys.argv[x] == "-inlinedvm":
        computeinlinedvm = bool(1);
    if sys.argv[x] == "-maxlen":
        _MAXLEN = int(sys.argv[x+1]);
    if sys.argv[x] == "-d":
        decompose = bool(1);
    if sys.argv[x] == "-b":
        buildgrammar = bool(1);
    if sys.argv[x] == "-actualavg":
        histmethod = ACTUALAVG;
    if sys.argv[x] == "-multi":
        histmethod = MULTI;
    if sys.argv[x] == "-randmulti":
        histmethod = RANDMULTI;
    if sys.argv[x] == "-dir":
        directory = sys.argv[x + 1];
        if directory[len(directory) - 1] != "\\" and directory[len(directory) - 1] != "/":
            print("error: terminate dir flag value with directory delimiter, i.e. end it in a slash (forward or backward depending on your system)");
            sys.exit(0);
    if sys.argv[x] == "-matrix":
        matrixfile = sys.argv[x + 1];
    if sys.argv[x] == "-paths":
        completedpathsfile = sys.argv[x + 1].split(',');
    if sys.argv[x] == "-i":
        if sys.argv[x + 1] == "s":
            stoichiometric = 1;
        if sys.argv[x + 1] == "a":
            stoichiometric = 0;    
    if sys.argv[x] == "-c":
        coefficients = bool(1);
    if sys.argv[x] == "-dmine":
        dmine = bool(1);
    if sys.argv[x] == "-log":
        loadlog = bool(1);
    if sys.argv[x] == "-e":
        _ADDENVIRON = bool(1);
    if sys.argv[x] == "-map":
        symbol_map_ = sys.argv[x + 1].split(',');    
    if sys.argv[x] == "-num":
        instancenum = int(sys.argv[x+1]);
    if sys.argv[x] == "-makeadj":
        makeadj = bool(1);
    if sys.argv[x] == "-sn":
        separateneuse = bool(1);
    if sys.argv[x] == "-cn":
        computeneuse = bool(1);
    if sys.argv[x] == "-an":
        aggregateneuse = bool(1);        
    if sys.argv[x] == "-bwan":
        buildadjneuse = bool(1);
    if sys.argv[x] == "-revtran": #[dir] [comma separated list of coefficients to use]
        if sys.argv[x + 1][len(sys.argv[x + 1]) - 1] != "\\" and sys.argv[x + 1][len(sys.argv[x + 1]) - 1] != "/":
            print("error: terminate dir command value with directory delimiter, i.e. end it in a slash (forward or backward depending on your system)");
            sys.exit(0);
            
        revtran_adj = LoadAndZeroAdj(sys.argv[x + 1], "matrix.txt");
        revtran_include = { };
        for item in sys.argv[x+2].split(','):
            revtran_include[int(item)] = 0;

        rev = ReverseTransform(sys.argv[x+1], revtran_adj, revtran_include);

        revfile = open(sys.argv[x+1] + "rev_matrix.txt", 'w');
        for row in rev:
            for col in row:
                revfile.write(str(col) + " ")
            revfile.write("\n");
        revfile.close();
        
        sys.exit();
    if sys.argv[x] == "-revtran2": #[dir] [comma separated list of coefficients to use] [in file name]

        if (os.path.exists(sys.argv[x+1] + "master.txt")):
            os.remove(sys.argv[x+1] + "master.txt");
        
        masterout = open(sys.argv[x+1] + "master.txt", 'w');
        masterin = open(sys.argv[x+1] + sys.argv[x+3], 'r');
                
        for line in masterin:
            if (os.path.exists(sys.argv[x+1] + "freqcount.txt")):
                os.remove(sys.argv[x+1] + "freqcount.txt");
            freqfile = open(sys.argv[x+1] + "freqcount.txt", 'w');
            cnt = 0;
            for item in line.split(','):
                if cnt >= 3:
                    freqfile.write(item + "\n");
                else:
                    if cnt != 2: masterout.write(item + " ");
                
                cnt += 1;
                
            freqfile.close();
            
            revtran_adj = LoadAndZeroAdj(sys.argv[x + 1], "matrix.txt");
            revtran_include = { };
            for item in sys.argv[x+2].split(','):
                revtran_include[int(item)] = 0;

            rev = ReverseTransform(sys.argv[x+1], revtran_adj, revtran_include);
            
            for row in rev:
                for col in row:
                    masterout.write(str(col) + " ")
                
            masterout.write("\n");

        masterin.close();
        masterout.close();
        
        sys.exit();

    
graph = None;
inttosymbol = None;
symboltoint = None;

if makeadj:
   graph, weightedgraph, inttosymbol, symboltoint = buildAdjacencyMatrixFromFile(directory, completedpathsfile[0], _ADDENVIRON, symbol_map_);    
         
if decompose:
    graph, flux_indexer, fluxmatrix, set_chains, set_cycles, out_chains, out_cycles, inttosymbol, symboltoint = Decompose(graph, 0, symbol_map_, inttosymbol, symboltoint);
    
if DVM:
    DVMSearchReal(set_chains, set_cycles, len(graph), out_chains, out_cycles);
    
if loadgraph:
    if stoichiometric:
        graph = LoadStoichiometricMatrix();
    else:
        graph = LoadAdjacencyMatrix(symbol_map_);

if loadfluxindexer:
    flux_indexer = LoadFluxIndexer();

if buildgrammar:
    GRAMMAR, _X, flux_indexer, reverseflux_indexer = BuildGrammar(graph,flux_indexer,0+1);    
    
if savegrammar:
    SaveGrammar(GRAMMAR, _X);

if loadgrammar:
    GRAMMAR, _X = LoadGrammar();

if coefficients:
    vectorlen = len(flux_indexer);
    histogram = [0.0 for i in range(vectorlen)];
    strcnt = 0;
    dvm = [];
    if loadlog:
        logfile = open(directory + "log.txt", 'r');
        cnt = 0;
        for line in logfile:
            lvals = line.strip("\r").strip("\n").split(',');
            if cnt == 0:
                strcnt = int(lvals[0]);
            elif cnt == 1:
                for i in range(vectorlen):
                    histogram[i] = float(lvals[i]);
            else:
                dvm.append([0 for i in range(vectorlen)]);
                for i in range(vectorlen):
                    dvm[len(dvm)-1][i] = int(lvals[i]);
            cnt += 1;
        logfile.close();

    pathdict = {};
    histogram, dvm, pathdict = ComputeCoefficientsAndDVM(flux_indexer,GRAMMAR,_X,graph,histogram,strcnt,dvm,completedpathsfile[0],symboltoint,None,None,None,_ADDENVIRON,pathdict,instancenum);
    
    csvheader = "";
    csvrec = "";
    
    file = open(directory + "freqcount.txt", 'w');
    
    #file.write(MakeString(histogram, 0, len(histogram)).strip(','));
    for i in range(0, len(histogram)):
        file.write(str(histogram[i]));
        file.write("\n");
        csvrec += (",", "")[len(csvrec) < 1] + str(histogram[i]);
        
    file.close();

    file = open(directory + "freqcount_flux.txt", 'w');
    
    for i in range(0, len(histogram)):
        file.write(reverseflux_indexer[i]);
        file.write("\n");
        csvheader += (",", "")[len(csvheader) < 1] + reverseflux_indexer[i].replace(",", " ");
        
    file.close();

    csvfile = open(directory + "out.csv", 'w');
    csvfile.write(csvheader + "\n");
    csvfile.write(csvrec + "\n");
    csvfile.close();
    
    file = open(directory + "dvm-inlinesearch.txt", 'w');
    for item in dvm: file.write(MakeString(item, 0, len(item)).strip(',') + "\n");
    file.close();
