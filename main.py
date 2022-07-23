import matplotlib.pyplot as plt
import math
import networkx as nx
from scipy import linalg as LA
from itertools import combinations
import collections
import random
from random import randint
import os, sys


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
    for i in range(len(comb_list)):
        z = []
        v = []
        w = []
        neighbours = []
        for j in range(k):
            # Find neighbourhood
            v = [comb_list[j]]
            neighbours = v + list(G.neighbors(comb_list[j]))  # neighborhood of v and v included
            z.append(neighbours)

        Int = list(set.intersection(*map(set, z)))

        H = G.subgraph(Int)
        cc = int(10 * math.sqrt(H.number_of_nodes())) + 1
        A = nx.adjacency_matrix(H)
        e_values, e_vecs = LA.eig(A.todense())
        sorted_indx = sorted(range(len(e_values)), key=lambda k: e_values[k], reverse=True)
        e_values_sorted = sorted(e_values, reverse=True)
        l_2 = e_values[sorted_indx[1]]  # 2nd largest eigenvalue

        index = sorted_indx[1]
        v_2 = e_vecs[sorted_indx[1]]  # neighborhood of v and v included
        v_2 = abs(v_2)
        sorted_indx_v2 = sorted(range(len(v_2)), key=lambda k: v_2[k], reverse=True)
        v_2 = sorted(v_2, reverse=True)
        sorted_vertices = []
        for j in range(len(sorted_indx)):
            sorted_vertices.insert(j, Int[sorted_indx_v2[j]])

        w = []
        w = sorted_vertices[0:c]
        minimum_neighbours = int(math.ceil(3 * c / 4))
        nodes = []
        common = []
        A = []
        for l in range(H.number_of_nodes()):
            A = H.neighbors(list(H.nodes())[l])
            A = list(A)
            common = list(set(A) & set(w))
            if (len(common)) >= minimum_neighbours:
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

    for i in range(len(comb_list)):
        z = []
        v = []
        w = []
        neighbours = []
        for j in range(k):
            # Find neighbourhood
            v = [comb_list[i][j]]
            neighbours = v + list(G.neighbors(comb_list[i][j]))  # neighborhood of v and v included
            z.append(neighbours)

        Int = list(set.intersection(*map(set, z)))
        H = G.subgraph(Int)
        cc = int(10 * math.sqrt(H.number_of_nodes())) + 1

        A = nx.adjacency_matrix(H)
        e_values, e_vecs = LA.eig(A.todense())
        sorted_indx = sorted(range(len(e_values)), key=lambda k: e_values[k], reverse=True)
        l_2 = e_values[sorted_indx[1]]  # 2nd largest eigenvalue
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
        minimum_neighbours = int(math.ceil(3 * c / 4))
        nodes = []

        common = []
        A = []
        for l in range(H.number_of_nodes()):
            A = H.neighbors(list(H.nodes())[l])
            A = list(A)
            common = list(set(A) & set(w))

            if (len(common)) >= minimum_neighbours:
                if not nodes:
                    nodes.append(list(H.nodes())[l])
                else:
                    count = 0
                    for nn in nodes:
                        if H.has_edge(nn, list(H.nodes())[l]) == True:
                            count = count + 1  # counting neighbours
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


def maximum_clique_set(G, k, file_clique, comb_list):
    M = []
    cliques = []
    L = []
    iter = 0
    z = []
    v = []
    w = []
    neighbours = []

    found = True
    finish = False
    flag_found = False
    first_time = 0
    failed_times = 0
    # Find neighbourhood
    if k == 1:
        flag_found = False
        v = int(''.join(str(i) for i in comb_list))
        a = list(G.neighbors(v))
    else:
        v = comb_list
        a = list(G.neighbors(comb_list))
    neighbours = [v] + a  # neighborhood of v and v included
    z.append(neighbours)
    Int = list(set.intersection(*map(set, z)))

    if not Int:
        r = 0
        return r

    H = G.subgraph(Int)
    c1 = 0
    for f in H.edges():
        c1 = c1 + 1
    non_edges_list = list(nx.non_edges(H))
    count = 0
    if len(non_edges_list) == 0:
        flag_found = True

    for j in H.edges():
        count = count + 1
    while len(non_edges_list) != 0:
        flag_found = False
        r_number1 = random.randint(0, len(non_edges_list) - 1)
        r_number2 = random.randint(0, 1)
        r_vertex = non_edges_list[r_number1][r_number2]
        if first_time == 0:
            H_delete = nx.Graph(H)
        else:
            H_delete = nx.Graph(H_delete)
        first_time = first_time + 1
        H_delete.remove_node(r_vertex)

        c1 = 0
        for f in H_delete.edges():
            c1 = c1 + 1  # count number of edges of H_delete
        V = H_delete.number_of_nodes()
        failed_times = failed_times + 1
        if c1 == V * (V - 1) / 2: #complete graph
            non_edges_list_after = list(nx.non_edges(H_delete))
            if found == False:
                fraction2 = abs((len(file_clique[1]) - failed_times) / len(file_clique[1]))
                finish = True
            break
        else: #not a complete graph
            non_edges_list = list(nx.non_edges(H_delete))
            found = False
            continue
    L.sort()
    Int.sort()
    L.append(Int)
    
    L = sorted(L, key=len, reverse=True)
    M = []

    if not L: #no cliques found
        M = []
        return M
    else:
        L_ = sorted(L, key=len, reverse=True)
        L_ = [item for L_ in L_ for item in L_]
        if finish:
            return [L_, fraction2]
        else:
            fraction1 = len(L_) / len(file_clique[1])
            return [L_, fraction1]


def maximum_clique_setk2(G, k, file_clique, comb_list):
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
    c1 = 0
    for f in H.edges():
        c1 = c1 + 1  # count number of edges of H_delete
    
    non_edges_list = list(nx.non_edges(H))

    failed_times = 0
    first_time = 0
    found = True
    flag_found = False
    finish = False
    sum_of_failed_fractions = 0
    if len(non_edges_list) == 0:
        flag_found = True
    while len(non_edges_list) != 0:
        flag_found = False
        r_number1 = random.randint(0,len(non_edges_list)-1)
        r_number2 = random.randint(0, 1)
        r_vertex = non_edges_list[r_number1][r_number2]
        if first_time == 0:
            H_delete = nx.Graph(H)
        else:
            H_delete = nx.Graph(H_delete)
        first_time = first_time + 1
        H_delete.remove_node(r_vertex)
        c1 = 0
        for f in H_delete.edges():
            c1 = c1 + 1  # count number of edges of H_delete
        V = H_delete.number_of_nodes()
        failed_times = failed_times + 1
        if c1 == V * (V - 1) / 2: #graph is complete
            if found == False:
                fraction2 = abs((len(file_clique[1]) - failed_times) / len(file_clique[1]))
                finish = True
            break
        else:
            non_edges_list = list(nx.non_edges(H_delete))
            found = False
            continue

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
        if finish:
            return [M_[0], fraction2]
        else:
            fraction1 = len(M_[0]) / len(file_clique[1])
            return [M_[0], fraction1]


def give_set_spectral(G, k, t, file_clique, RANDOM_COMBINATIONS, sum_of_failed_vertices):
    found_clique = 0
    counter_failed = 0
    sum_of_failed_fraction = 0
    list_of_failed_vertices = []

    if k == 1:
        for i in range(len(file_clique[1])):
            largest_clique = spectral_max_clique_set(G, k, t, [file_clique[1][i]])

            if sorted(largest_clique) == sorted(file_clique[1]):
                found_clique = 1
                fraction = len(largest_clique) / len(file_clique[1])
                return [found_clique, fraction]
            elif len(largest_clique) != len(file_clique[1]):
                counter_failed += 1
                failed_vertices = (len(file_clique[1]) - len(largest_clique))
                list_of_failed_vertices.append(failed_vertices)
                fraction = len(largest_clique) / len(file_clique[1])
                sum_of_failed_fraction += fraction
                continue
            else:
                counter_failed += 1
                continue
        if found_clique == 0:
            found_clique = 0
            return [found_clique, sum_of_failed_vertices, counter_failed, fraction]
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
                    if len(largest_clique) == len(file_clique[1]): #max clique found
                        found_clique = 1
                        fraction = len(largest_clique) / len(file_clique[1])
                        return [found_clique, fraction]
                    else: #failed to find max clique
                        failed_vertices = 0
                        failed_vertices = (len(file_clique[1]) - len(largest_clique))
                        list_of_failed_vertices.append(failed_vertices)
                        fraction = len(largest_clique) / len(file_clique[1])
                        sum_of_failed_fraction += fraction
                        continue

                if found_clique == 0:
                    found_clique = 0
                    return [found_clique, sum_of_failed_fraction, counter_failed, fraction]
            else:
                counter = counter + 1

def is_subclique(G,nodelist):
    H = G.subgraph(nodelist)
    n = len(nodelist)
    return H.size() == n*(n-1)/2

def give_set_max_clique(G, p, k, file_clique, fractions_list, failed_sum_list):
    r_comb = 20
    if k == 1:
        failed = False
        counter_failed = 0
        found_clique = 0
        total_combinations = math.comb(len(file_clique[1]), k)
        x = [randint(1, total_combinations) for p in range(0, r_comb)]
        x = sorted(x)
        counter2 = 0
        for i in range(len(file_clique[1])):
            counter2 = counter2 + 1
            largest_clique = maximum_clique_set(G, k, file_clique, [file_clique[1][i]])
            fraction = largest_clique[1]

            if len(largest_clique[0]) == len(file_clique[1]):

                if failed:
                    avg_fraction = failed_sum_list[-1] / counter_failed
                    fractions_list.append(avg_fraction)
                    return [found_clique, fractions_list]
                found_clique = 1
                fractions_list.append(fraction)
                return [found_clique, fractions_list]
            else:
                failed = True
                counter_failed = counter_failed + 1
                failed_sum_list.append(fraction)
                continue

        if found_clique == 0:
            avg_fraction = sum(failed_sum_list) / counter_failed
            fractions_list.append(avg_fraction)
           return [found_clique, fractions_list]
    else:
        total_combinations = math.comb(len(file_clique[1]), k)
        x = [randint(1, total_combinations) for p in range(0, r_comb)]
        x = sorted(x)
        counter = 1
        x_i = 0
        found_clique = 0
        counter_failed = 0
        for comb_list in combinations(file_clique[1], k):
            if x_i > r_comb:
                break
            if counter == x[x_i]:
                x_i += 1
                for i in comb_list:
                    largest_clique = maximum_clique_setk2(G, k, file_clique, comb_list)
                    fraction = largest_clique[1]
                    if len(largest_clique[0]) == len(file_clique[1]):
                        found_clique = 1
                        fractions_list.append(fraction)
                        return [found_clique, fractions_list]
                    else:
                        counter_failed = counter_failed + 1
                        failed_sum_list.append(fraction)
                        continue;
                if found_clique == 0:
                    avg_fraction = sum(failed_sum_list) / counter_failed
                    fractions_list.append(avg_fraction)
                    return [found_clique, fractions_list]
            else:
                counter = counter + 1
                
def plot_experiments_failure(n, a, k, count_spectral_false, count_max_clique_false, probability_vector):
    plt.close()
    plt.xlabel('p')
    plt.ylabel('Failure Probability')
    if a == 1/3:
        plt.title('n={}, a=1/3, k={}'.format(n, k))
    elif a == 2/3:
        plt.title('n={}, a=2/3, k={}'.format(n, k))
    else:
        plt.title('n={}, a=1, k={}'.format(n, k))

    plt.plot(probability_vector, count_spectral_false, "-r", label="Spectral-Max-Clique")
    plt.plot(probability_vector, count_max_clique_false, "-b", label="Maximum-Clique")
    plt.grid(True, color="grey", linewidth="0.5", linestyle="-.", label=True)
    plt.legend()
    plt.savefig("a{}k{}.png".format("%.2f" % a, k), dpi=300)

                
def plot_experiments_fractions(n, a, k, fraction_spectral, fraction_max_clique, probability_vector):

    plt.close()
    plt.xlabel('p')
    plt.ylabel('Approximation guarantee (ag)')
    if a == 1/3:
        plt.title('n={}, a=1/3, k={}'.format(n, k))
    elif a == 2/3:
        plt.title('n={}, a=2/3, k={}'.format(n, k))
    else:
        plt.title('n={}, a=1, k={}'.format(n, k))

    plt.plot(probability_vector, fraction_spectral, "-r", label="Spectral-Max-Clique")
    plt.plot(probability_vector, fraction_max_clique, "-b", label="Maximum-Clique")
    plt.grid(True, color="grey", linewidth="0.5", linestyle="-.", label=True)
    plt.legend()
    plt.savefig("fraction_a{}k{}.png".format("%.2f" % a, k), dpi=300)


def run_tests(n, a, k, probabilities):
    sum_all_list = []
    sum_all_list_max_clique = []
    max_fraction_list = []
    fractions_list = []
    m = math.floor(n ** a)
    experiment_count = 2 #times each experiment will be executed (e.g 100 to get a percentage)

    percentage_list_spectral = []
    percentage_list_max_clique = []
    RANDOM_COMBINATIONS = 20

    spectral_max_file = open("experiments_Spectral_Max_Clique.txt", "a")
    spectral_max_file.write(
        "\nn=" + str(n) + ", a=" + str(a) + ", k=" + str(k) + ", Probability Vector =" + str(probabilities) + "\n")
    spectral_max_file.write("\n")

    maximum_clique_file = open("experiments_Maximum_Clique.txt", "a")
    maximum_clique_file.write(
        "n=" + str(n) + ", a=" + str(a) + ", k=" + str(k) + ", Probability Vector =" + str(probabilities) + "\n")
    maximum_clique_file.write("\n")

    for prob in probabilities:
        false_prob_spectral = []
        true_spectral = 0
        false_spectral = 0

        flag_max_clique = False
        flag_spectral = False
        false_prob_max_clique = []
        true_max_clique = 0
        false_max_clique = 0

        sum_of_failed_vertices = 0
        counter_succeed = 0
        fractions_all_list = []
        fractions_list = []
        for i in range(0, experiment_count):
            failed_sum_list = []
            pos = [prob] * m
            t = int(n * prob)
            mlabels, G = nx.general_random_intersection_graph(n, m, pos)
            file_clique = printToFile(mlabels)

            found_or_not = give_set_spectral(G, k, t, file_clique, RANDOM_COMBINATIONS,sum_of_failed_vertices)

            if len(found_or_not) == 4: #meaning algorithm did not found the max clique
                fractions_all_list.append(found_or_not[3])
            else:
                counter_succeed += 1
            if found_or_not == 1:
                true_spectral += 1
                flag_spectral = True
            r = give_set_max_clique(G, prob, k, file_clique, fractions_list, failed_sum_list)

            if r[0] == 1:failed_sum_list
                true_max_clique += 1
                flag_max_clique = True
            else:
                flag_max_clique = False
                false_max_clique += 1

            false_prob_max_clique.append(flag_max_clique)
            false_prob_spectral.append(flag_spectral)

        max_fraction = 0
        max_fraction_list.append(max_fraction)
        sum_of_failed_fractions = sum(fractions_all_list)
        sum_of_succeded = counter_succeed
        sum_all = (sum_of_failed_fractions + sum_of_succeded) / experiment_count
        sum_all = round(sum_all, 2)
        sum_max_clique_fraction = round(sum(r[1]) / experiment_count, 2)
        spectral_max_file.write("p= " + str(prob) + " Failure Percentage: " + str(false_spectral) + "\n")
        maximum_clique_file.write("p= " + str(prob) + "Failure Percentage: " + str(false_max_clique) + "\n")
        spectral_max_file.write("\n")
        maximum_clique_file.write("\n")

        percentage_list_spectral.append(false_spectral)
        percentage_list_max_clique.append(false_max_clique)
        sum_all_list.append(sum_all)
        sum_all_list_max_clique.append(sum_max_clique_fraction)

    spectral_max_file.write("k=" + str(k) + ", a=" + str(a) + ", p=" + str(probabilities))
    spectral_max_file.write("\n")
    spectral_max_file.write("Sum spectral fraction: " + str(sum_all_list))
    spectral_max_file.write("\n")
    spectral_max_file.write("\n")

    maximum_clique_file.write("k=" + str(k) + ", a=" + str(a) + ", p=" + str(probabilities))
    maximum_clique_file.write("\n")
    maximum_clique_file.write("Sum max clique fraction: " + str(sum_all_list_max_clique))
    maximum_clique_file.write("\n")
    plot_experiments_failure(n, a, k, sum_all_list_max_clique, sum_all_list_max_clique, probabilities)
    plot_experiments_fractions(n, a, k, percentage_list_spectral, percentage_list_max_clique, probabilities)
