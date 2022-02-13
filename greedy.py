def greedy_clique(G):
    start = timeit.default_timer()
    degrees = sorted(G.degree, key=lambda x: x[1], reverse=True)
    M = []

    for i in range(len(degrees)):
        if not M:  # if M is empty
            M.append(degrees[i][0])
            continue
        c = 0
        for k in M:
            if G.has_edge(degrees[i][0], k):
                c = c + 1
        if c == len(M):
            M.append(degrees[i][0])

    return M

def mono_clique(G):
    intersections = []
    for edge in G.edges():
        neighbours = []
        z = []
        neighbours.append([edge[0]] + list(G.neighbors(edge[0])))
        neighbours.append([edge[1]] + list(G.neighbors(edge[1])))
        intersections.append(list(set.intersection(*map(set, neighbours))))

    D = sorted(intersections, key=len, reverse=True)

    for i in D:
        count = 0
        H = G.subgraph(i)
        for e in H.edges():
            count = count + 1

        V = H.number_of_nodes()
        if count == V * (V - 1) / 2:
            return i
        else:
            continue;
