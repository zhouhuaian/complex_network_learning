# 计算度分布
def get_pdf(G):
    all_k = [G.degree[i] for i in G.nodes]  # 获得每个节点的度
    k = list(set(all_k))  # 去掉重复的度值
    N = len(G.nodes)  # 节点数目

    pk = []  # 度分布
    sorted_k = sorted(k)
    # 对于每个度值，统计节点所占比例，即该度值的概率
    for ki in sorted_k:
        cnt = 0
        for i in G.nodes:
            if G.degree[i] == ki:
                cnt += 1
        pk.append(cnt / N)
    
    return sorted_k, pk