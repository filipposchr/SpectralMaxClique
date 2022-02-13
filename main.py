import matplotlib.pyplot as plt
import math
import networkx as nx
from scipy import linalg as LA
from itertools import combinations
import timeit
import collections
import random
from random import randint

def printToFile(mlabels):
    reverse_list = map(lambda x: x[::-1], mlabels)
    result = {}

    for group in reverse_list:
        for item in group[1:]:
            result.setdefault(item, []).append(group[0])

    result = collections.OrderedDict(sorted(result.items()))
    count = {}
    lr = []

    max = 0
    for key, values in result.items():
        count[key] = len(values)
        if len(values) > max:
            max = len(values)
            max_clique = values
        lr.append(len(values))
    lr = sorted(lr, reverse=True)

    with open('labelsNodes.txt', 'w') as f:
        for key, value in result.items():
            f.write('Label {0},  Nodes : {1}\n'.format(key, value))

    return lr[0], max_clique

def spectral_max_clique_set(G, k, c, comb_list):
    M = []
    cliques = []
    geit = 0
    for i in range(len(comb_list)):
        z = []
        v = []
        w = []
        neighbours = []
        for j in range(k):
            # Find neighbourhood
            v = [comb_list[j]]
            neighbours = v + list(G.neighbors(comb_list[j]))#neighborhood of v and v included
            z.append(neighbours)

        Int = list(set.intersection(*map(set, z)))

        H = G.subgraph(Int)
        cc = int(10 * math.sqrt(H.number_of_nodes())) + 1
        A = nx.adjacency_matrix(H)
        e_values, e_vecs = LA.eig(A.todense())
        sorted_indx = sorted(range(len(e_values)), key=lambda k: e_values[k], reverse=True)
        e_values_sorted = sorted(e_values, reverse=True)
        l_2 = e_values[sorted_indx[1]] #2nd largest eigenvalue
        korifi_2i = Int[sorted_indx[1]]

        index = sorted_indx[1]
        v_2 = e_vecs[sorted_indx[1]]#neighborhood of v and v included
        v_2 = abs(v_2)
        sorted_indx_v2 = sorted(range(len(v_2)), key=lambda k: v_2[k], reverse=True)
        v_2 = sorted(v_2, reverse=True)
        sorted_vertices = []
        for j in range(len(sorted_indx)):
            sorted_vertices.insert(j, Int[sorted_indx_v2[j]])

        w = []
        w = sorted_vertices[0:c]
        toul_geit = int(math.ceil(3 * c / 4))
        nodes = []
        koina = []
        A = []
        for l in range(H.number_of_nodes()):
            A = H.neighbors(list(H.nodes())[l])
            A = list(A)
            koina = list(set(A) & set(w))
            if (len(koina)) >= toul_geit:
                if not nodes:
                    nodes.append(list(H.nodes())[l])
                else:
                    count = 0
                    for nn in nodes:
                        if (H.has_edge(nn, list(H.nodes())[l]) == True):
                            count = count + 1
                    if count >= (len(nodes)):
                        nodes.append(list(H.nodes())[l])
        if not nodes:
            continue
        else:
            cliques.append(nodes)
    cliques = sorted(cliques, key=len, reverse=True)
    try:
        largest_clique = sorted(cliques[0])
    except IndexError:
        largest_clique = []
        return largest_clique
    else:
        return largest_clique

def spectral_max_clique_set_k2(G, k, c, comb_list):
    M = []
    cliques = []
    geit = 0
    for i in range(len(comb_list)):
        z = []
        v = []
        w = []
        neighbours = []
        for j in range(k):
            # Find neighbourhood
            v = [comb_list[i][j]]
            neighbours = v + list(G.neighbors(comb_list[i][j])) #neighborhood of v and v included
            z.append(neighbours)

        Int = list(set.intersection(*map(set, z)))
        H = G.subgraph(Int)
        cc = int(10 * math.sqrt(H.number_of_nodes())) + 1

        A = nx.adjacency_matrix(H)
        e_values, e_vecs = LA.eig(A.todense())
        sorted_indx = sorted(range(len(e_values)), key=lambda k: e_values[k], reverse=True)
        l_2 = e_values[sorted_indx[1]]  #2nd largest eigenvalue
        korifi_2i = Int[sorted_indx[1]]

        index = sorted_indx[1]
        v_2 = e_vecs[sorted_indx[1]]

        v_2 = abs(v_2)
        sorted_indx_v2 = sorted(range(len(v_2)), key=lambda k: v_2[k], reverse=True)
        v_2 = sorted(v_2, reverse=True)

        sorted_vertices = []
        for j in range(len(sorted_indx)):
            sorted_vertices.insert(j, Int[sorted_indx_v2[j]])
        w = []
        w = sorted_vertices[0:c]
        toul_geit = int(math.ceil(3 * c / 4))
        nodes = []

        koina = []
        A = []
        for l in range(H.number_of_nodes()):
            A = H.neighbors(list(H.nodes())[l])
            A = list(A)
            koina = list(set(A) & set(w))

            if (len(koina)) >= toul_geit:
                if not nodes:
                    nodes.append(list(H.nodes())[l])
                else:
                    count = 0
                    for nn in nodes:
                        if (H.has_edge(nn, list(H.nodes())[l]) == True):
                            count = count + 1 #counting neighbours
                    if count >= (len(nodes)):
                        nodes.append(list(H.nodes())[l])

        if not nodes:
            continue
        else:
            cliques.append(nodes)

    cliques = sorted(cliques, key=len, reverse=True)
    found = 0
    try:
        largest_clique = sorted(cliques[0])
    except IndexError:
        largest_clique = []
        return largest_clique
    else:
        return largest_clique

def maximum_clique_set(G, k, comb_list):
    M = []
    cliques = []
    geit = 0
    L = []
    ll = []
    iter = 0
    z = []
    v = []
    w = []
    neighbours = []

    # Find neighbourhood
    if k == 1:
        v = int(''.join(str(i) for i in comb_list))
        a = list(G.neighbors(v))
    else:
        v = comb_list
        a = list(G.neighbors(comb_list))
    neighbours = [v] + a  #neighborhood of v and v included
    z.append(neighbours)

    for prosk_v in G.neighbors(v):
        geitonikes_v = []
        geitonikes_v = [prosk_v] + list(G.neighbors(prosk_v))

        if set(geitonikes_v) == set(neighbours):
            ll.append(prosk_v)
    Int = list(set.intersection(*map(set, z)))

    if not Int:
        r = 0
        return r

    H = G.subgraph(Int)
    count = 0

    for j in H.edges():
        count = count + 1

    V = H.number_of_nodes()
    if count == V * (V - 1) / 2:
        L.sort()
        Int.sort()

        if Int not in L:
            L.append(Int)
        else:
            r = 0
            return r
    iter = iter + 1
    Lnew = []

    L = sorted(L, key=len, reverse=True)
    Y = []
    M = []
    count = 0
    
    if k > 1:
        for i in range(len(L)):
            perm = []
            perm = combinations(list(L[i]), 2)
            perm = list(perm)

            for combination in perm:
                a = sorted(combination, key=lambda x: int(x), reverse=True)
                b = sorted(combination, key=lambda x: int(x), reverse=False)

                if ((tuple(combination) not in Y) or (tuple(b) not in Y)):
                    Y.extend(perm)
                    M.append(L[i])
                    break
        if not M:
            print("No Cliques Found.")
        else:
            count = 0
            for i in range(len(M)):
                for x in M[i]:
                    count = count + 1
            M_ = sorted(M, key=len, reverse=True)
            return M_[0]
    else:
        if not L:
            print("No Cliques Found.")
            L = []
            return L
        else:
            L_ = sorted(L, key=len, reverse=True)
            L_ = [item for L_ in L_ for item in L_]
            return L_


def maximum_clique_setk2(G, k, comb_list):
    L = []
    z = []

    for j in range(k):
        v = [comb_list[j]]
        neighbours = v + list(G.neighbors(comb_list[j]))  # neighborhood of v and v included
        z.append(neighbours)

    Int = list(set.intersection(*map(set, z)))

    if not Int:
        M = []
        return M

    H = G.subgraph(Int)

    count = 0
    for f in H.edges():
        count = count + 1

    V = H.number_of_nodes()
    if count == V * (V - 1) / 2:
        Int.sort()
        L.append(Int)

    L = sorted(L, key=len, reverse=True)
    Y = []
    M = []
    if not L:
        M = []
        return M

    for i in range(len(L)):
        perm = []
        perm = combinations(list(L[i]), 2)
        perm = list(perm)

        for combination in perm:
            a = sorted(combination, key=lambda x: int(x), reverse=True)
            b = sorted(combination, key=lambda x: int(x), reverse=False)

            if ((tuple(combination) not in Y) or (tuple(b) not in Y)):
                Y.extend(perm)
                M.append(L[i])
                break
    if not M:
        M = []
        return M
    else:
        M_ = sorted(M, key=len, reverse=True)
        return M_[0]


def give_set(G, k, t, file_clique, RANDOM_COMBINATIONS):
    found_clique = 0
    if k == 1:
        for i in range(len(file_clique[1])):
            largest_clique = spectral_max_clique_set(G, k, t, [file_clique[1][i]])

            if sorted(largest_clique) == sorted(file_clique[1]):
                found_clique = 1
                return found_clique
            else:
                continue;
        if found_clique == 0:
            return found_clique
    else:
        total_combinations = math.comb(len(file_clique[1]), k)
        x = [randint(1, total_combinations) for p in range(0, RANDOM_COMBINATIONS)]
        x = sorted(x)
        step = int(math.ceil(total_combinations / RANDOM_COMBINATIONS))
        counter = 1
        x_i = 0
        for comb_list in combinations(file_clique[1], k):
            if x_i > RANDOM_COMBINATIONS:
                break
            if (counter == x[x_i]):
                x_i += 1
                for i in comb_list:
                    largest_clique = spectral_max_clique_set_k2(G, k, t, [comb_list])
                    if sorted(largest_clique) == sorted(file_clique[1]):
                        found_clique = 1
                        return found_clique
                    else:
                        continue;
                if found_clique == 0:
                    found_clique = 0
                    return found_clique
            else:
                counter = counter + 1

def give_set_featureFind(G,p, k, file_clique):
    found_clique = 0
    RANDOM_COMBINATIONS = 20
    if k == 1:
        total_combinations = math.comb(len(file_clique[1]), k)
        x = [randint(1, total_combinations) for p in range(0, RANDOM_COMBINATIONS)]
        x = sorted(x)
        counter = 1
        x_i = 0
        for i in range(len(file_clique[1])):
            largest_clique = maximum_clique_set(G, k, [file_clique[1][i]])
            if sorted(largest_clique) == sorted(file_clique[1]):
                found_clique = 1
                return found_clique
            else:
                continue;
        if found_clique == 0:
            return found_clique
    else:
        total_combinations = math.comb(len(file_clique[1]), k)
        x = [randint(1, total_combinations) for p in range(0, RANDOM_COMBINATIONS)]
        x = sorted(x)
        counter = 1
        x_i = 0
        for comb_list in combinations(file_clique[1], k):
            if x_i > RANDOM_COMBINATIONS:
                break
            if (counter == x[x_i]):
                x_i += 1
                for i in comb_list:
                    largest_clique = maximum_clique_setk2(G, k, comb_list)
                    if sorted(largest_clique) == sorted(file_clique[1]):
                        found_clique = 1
                        return found_clique
                    else:
                        continue;
                if found_clique == 0:
                    return found_clique
            else:
                counter = counter + 1


def run_spectral_tests():
    RANDOM_COMBINATIONS = 20
    text_file = open("expirements_Spectral_Max_Clique.txt", "a")
    text_file.write("\nn=" + str(n) + ", a=" + str(a) + ", k=" + str(k) + ", pinakas_pith =" + str(pinakas_pith2) + "\n")
    for pith in pinakas_pith:
        pith_lathous_alon = []
        swsta_alon = 0

        for i in range(0,20):
            pos = [pith] * m
            pithan.append(pith)
            t = int(n * pith)
            mlabels, G = nx.general_random_intersection_graph(n, m, pos)
            file_clique = printToFile(mlabels)
            found_or_not = give_set(G, k, t, file_clique, RANDOM_COMBINATIONS)

            if found_or_not == 1:
                swsta_alon = swsta_alon + 1
                flag_alon = True
            elif found_or_not == 0:
                flag_alon = False

            pith_lathous_alon.append(flag_alon)
        text_file.write("p= " + str(pith) + " pith_lathous: " + str(pith_lathous_alon) + "\n")
        text_file.write("swsta alon : " + str(swsta_alon) + "\n")

    text_file.write("\n")

def run_featurefind_tests():
    RANDOM_COMBINATIONS = 20
    text_file = open("expirements_Maximum_Clique.txt", "a")
    text_file.write("n=" + str(n) + ", a=" + str(a) + ", k=" + str(k) + ", pinakas_pith =" + str(pinakas_pith2) + "\n")
    for pith in pinakas_pith:
        flag_featurefind = False
        pith_lathous_featurefind = []
        swsta_featurefind = 0

        for i in range(0, 20):
            pos = [pith] * m
            pithan.append(pith)
            t = int(n * pith)

            mlabels, G = nx.general_random_intersection_graph(n, m, pos)
            file_clique = printToFile(mlabels)
            if k == 1:
                r = give_set_featureFind(G, pith, k, file_clique)
            else:
                r = give_set_featureFind(G, pith, k, file_clique)

            if k == 1:
                print("File clique: ", file_clique[1])
                print("r: " ,r)
                print("File clique: ", sorted(file_clique[1]))
                if r == 1:
                    swsta_featurefind = swsta_featurefind + 1
                    flag_featurefind = True
                else:
                    flag_featurefind = False
            else:
                if r == 1:
                    swsta_featurefind = swsta_featurefind + 1
                    flag_featurefind = True
                else:
                    flag_featurefind = False
            print("flag: ", flag_featurefind)
            pith_lathous_featurefind.append(flag_featurefind)

        text_file.write("p= " + str(pith) + " pith_lathous: " + str(pith_lathous_featurefind) + "\n")
        text_file.write("swsta ff : " + str(swsta_featurefind) + "\n")
    text_file.write("\n")

"""
"""

n = 1000
a = 1/3
k = 5
m = math.floor(n ** a)


pinakas_pith = [0.23, 0.25, 0.28, 0.3]
pinakas_pith2 =[0.23, 0.25, 0.28, 0.3]
pithan = []

run_spectral_tests()
run_featurefind_tests()

"""
"""
"""
pos = [pinakas_pith[0]] * m
mlabels, G = nx.general_random_intersection_graph(n, m, pos )
c = int(n * pinakas_pith[0])

spectral_max_clique(G,k,c)
larg_cliq = printToFile(mlabels)
print(larg_cliq)

#show_n_m()
#plt.subplot(111)
#plt.text(0, 0, 'n= ' + str(n) + ',m= ' + str(m) + ',a= ' + str(a  + ",p= " + str(p), fontsize=13, color='red',  transform=plt.gcf().transFigure)
#plt.show() 

"""
"""
pith = [0.008, 0.009, 0.01, 0.011, 0.0115, 0.012, 0.0125, 0.013, 0.0135, 0.014, 0.015, 0.016, 0.0165, 0.017, 0.0175,
        0.018, 0.019, 0.02]
ff = [0, 6, 4, 18, 27, 40, 62, 76, 74, 92, 100, 100, 100, 100, 100, 100, 100, 100]
alon = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 21, 43, 62, 84, 86, 100, 100]

plt.xlabel('p')
plt.ylabel('Failure Probability')
plt.title('n=3000, a=1, k=8')
plt.plot(pith, alon , "-r", label="Spectral-Max-Clique")
plt.plot(pith, ff, "-b", label="Maximum-Clique")
plt.grid(True, color = "grey", linewidth = "0.5", linestyle = "-.", label=True)
plt.legend()
plt.savefig("a1k8.png", dpi=300)
"""
"""
# print("Length Largest Clique Alon found: ", len_clique_alon)
# counter = 1


for i in range(0,100):
    flag = True
    mlabels, G = nx.general_random_intersection_graph(n, m, pos)
    ss1 = timeit.default_timer()
    #a = spectral_max_clique(G,1,3* int(avg_clique_size))
    #a = greedy_clique(G)
    #a = featureFind3(G,n,p,k)
    a = mono_clique(G)
    ss2 = timeit.default_timer()
    b = printToFile(mlabels)
    rr = ss2 - ss1
    #print(rr)
    time = time + rr
    if a != b:
        times.append(0)
        flag = False
    else:
        times.append(1)
    print(flag)
    #print(counter)
    #counter = counter + 1

#print("counter :" , counter)
#print("Xronos ektelesis: ", time / 100)
#print("Ta 0 einai: ", times.count(0))

#larg_cliq = printToFile(mlabels)
#avg_clique_size = int(n * p)
#spectral_max_clique(G,2,avg_clique_size)
#marked_featureFind(G,n,p,1)
#greedy_clique(G)
#mono_clique(G)
#greedy_clique_opt(G)

#Add planted clique
#clique_nodes = list(range(200, 227))
#combina = list(combinations(clique_nodes, 2))
#print("Combina: ", list(combina))
#G.add_edges_from(list(combina))

#print("(Node, Label): ", mlabels)

plt.subplot(111)

#show_n_only_with_edges()
#show_n_m()
#show_n()

#plt.savefig("1.png")
plt.text(0, 0, 'n= ' + str(n) + ',m= ' + str(m) + ',a= ' + str(a) + ",p= " + str(p), fontsize=13, color='red',  transform=plt.gcf().transFigure)
plt.show()

pith_lathous = list()
pithan = list()
a = 1
times = []
n = 100
m = math.floor(n ** a)
e = 0.01
p = 0.02

for o in range(0,1):
    pos = [p] * m
    mlabels, G = nx.general_random_intersection_graph(n, m, pos)
    times = []
    print("P: ", p)
    pithan.append(p)

    for i in range(0,10):

        flag = True
        mlabels, G = nx.general_random_intersection_graph(n, m, pos)
        ss1 = timeit.default_timer()
        #a = spectral_max_clique(G,1,3* int(avg_clique_size))
        #a = greedy_clique(G)
        #a = featureFind3(G,n,p,k)
        a = mono_clique(G)
        ss2 = timeit.default_timer()
        b = printToFile(mlabels)

        if a != b:
            times.append(0)
            flag = False
        else:
            times.append(1)
        print(flag)

    print("**********************", times.count(0))
    pith_lathous.append(times.count(0))
    p = p + e


plt.xlabel('p')
plt.ylabel('%Λάθους')
plt.title('Maximum-Clique')
plt.plot(pithan, pith_lathous, "-r", label="α = 1/3")
"""
