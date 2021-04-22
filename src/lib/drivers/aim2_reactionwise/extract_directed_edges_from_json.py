from networkx.readwrite import json_graph
import json
path = '/home/burkhart/Software/reticula/data/aim2/input/reactome_pathway_hierarchy.json'
f = open(path)
data = json.load(f)
f.close()

print(data[0])
print(data[0]['stId'])
print(data[0]['children'])
children = data[0]['children']
for i in children:
    print(i['stId'])
    print(i['children'])

