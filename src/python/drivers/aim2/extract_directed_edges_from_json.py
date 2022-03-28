from networkx.readwrite import json_graph
import json

#The .json dumped from Reactome is malformed and cannot be consumed in its original state.
#The below stack trace shows the library errors-out on raw file. Add a parent attribute
#such as 'Chemistry of Life' to put the text in a compliant/consumable state.
#Traceback (most recent call last):
#File "<input>", line 8, in <module>
#File "/home/burkhart/miniconda3/envs/reticula/python/python3.7/site-packages/networkx/readwrite/json_graph/tree.py", line 141, in tree_graph
#root = data[id_]
#TypeError: list indices must be integers or slices, not str

path = '/home/burkhart/Software/reticula/data/aim2/input/reactome_pathway_hierarchy.json'
f = open(path,'r')
data = json.load(f)
f.close()
H = json_graph.tree_graph(data, attrs=dict(id='stId', children='children'))

path = '/home/burkhart/Software/reticula/data/aim2/input/pathway_reaction_id_edges.txt'
f = open(path,'w')
for e in H.edges:
    parent_id = e[0]
    child_id = e[1]
    f.write("{child_id}\t{parent_id}\n".format(child_id=child_id, parent_id=parent_id))
f.close()



