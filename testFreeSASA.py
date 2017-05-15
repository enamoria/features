import freesasa

structure = freesasa.Structure("1ubq.pdb")
result = freesasa.calc(structure)
area_classes = freesasa.classifyResults(result, structure)

print "Total : %.2f A2" % result.totalArea()
for key in area_classes:
    print key, ": %.2f A2" % area_classes[key]