import features as F
from sequence import sequence

# Method getASAList moved to archived.py
# Using ASA method from features

def generateInputFromELM(elm):
    x = 0

def extractAll():
    f_out = open("features_output", "w")

    f = open("phosphoELM_all_2015-04.dump", "r")
    filestr = f.readline()

    while True:
        strTemp = f.readline()
        if strTemp == "":
            break

        strList = strTemp.split("\t")
        strList = strList[0:3]

        corpus = sequence(strList[1], strList[0], strList[2], window)

        f_out.write(str(F.feature_1_shannon_entropy(corpus)))
        f_out.write(',')
        f_out.write(str(F.feature_2_relative_entropy(corpus)))
        f_out.write(',')
        f_out.write(str(F.feature_3_information_gain(corpus)))
        f_out.write(',')

        ###### ASA ######
        f_ASA = open("./RVP_seq_output/" + strList[0] + ".out", "r")
        # ASAList = getASAList(f_ASA, strList[2])
        ASAList = F.feature_4_12_ASA(corpus)
        f_ASA.close()

        f_out.write(",".join(str(x) for x in ASAList))

        f_out.write(str(F.feature_13_102_overlapping_properties(corpus)))
        f_out.write(',')

        f_out.write(",".join(str(x) for x in F.feature_103_106_hydrophobicity(corpus)))
        f_out.write(',')

        f_out.write(",".join(str(x) for x in F.feature_287_433_CTD(corpus)))
        f_out.write(',')

        f_out.write(",".join(str(x) for x in F.feature_434_493_socn(corpus)))
        f_out.write(',')

        f_out.write(",".join(str(x) for x in F.feature_494_593_quasi(corpus)))

        f_out.write('\n')

def extractPeptide(peptide, accession, protein, phosPos, window = 9):
    f_out = open(peptide + ".out", "w")

    # f = open("phosphoELM_all_2015-04.dump", "r")
    # filestr = f.readline()

    # while True:
    #     strTemp = f.readline()
    #     if strTemp == "":
    #         break
    #
    #     strList = strTemp.split("\t")
    #     strList = strList[0:3]
    #
    #     corpus = sequence(strList[1], strList[0], strList[2], window)

    corpus = sequence(protein, accession, phosPos, window)

    f_out.write(str(F.feature_1_shannon_entropy(corpus)))
    f_out.write(',')

    f_out.write(str(F.feature_2_relative_entropy(corpus)))
    f_out.write(',')

    f_out.write(str(F.feature_3_information_gain(corpus)))
    f_out.write(',')

    ###### ASA ######
    f_ASA = open("./RVP_seq_output/" + accession + ".out", "r")
    ASAList = F.feature_4_12_ASA(corpus)
    f_ASA.close()

    f_out.write(",".join(str(x) for x in ASAList))
    f_out.write(',')

    f_out.write(str(F.feature_13_102_overlapping_properties(corpus)))
    f_out.write(',')

    f_out.write(",".join(str(x) for x in F.feature_103_106_hydrophobicity(corpus)))
    f_out.write(',')

    f_out.write(",".join(str(x) for x in F.feature_287_433_CTD(corpus)))
    f_out.write(',')

    f_out.write(",".join(str(x) for x in F.feature_434_493_socn(corpus)))
    f_out.write(',')

    f_out.write(",".join(str(x) for x in F.feature_494_593_quasi(corpus)))

string = "MAEMGSKGVTAGKIASNVQKKLTRAQEKVLQKLGKADETKDEQFEQCVQNFNKQLTEGTRLQKDLRTYLASVKAMHEASKKLSECLQEVYEPEWPGRDEANKIAENNDLLWMDYHQKLVDQALLTMDTYLGQFPDIKSRIAKRGRKLVDYDSARHHYESLQTAKKKDEAKIAKPVSLLEKAAPQWCQGKLQAHLVAQTNLLRNQAEEELIKAQKVFEEMNVDLQEELPSLWNSRVGFYVNTFQSIAGLEENFHKEMSKLNQNLNDVLVSLEKQHGSNTFTVKAQPSDNAPEKGNKSPSPPPDGSPAATPEIRVNHEPEPASGASPGATIPKSPSQLRKGPPVPPPPKHTPSKEMKQEQILSLFDDAFVPEISVTTPSQFEAPGPFSEQASLLDLDFEPLPPVASPVKAPTPSGQSIPWDLWEPTESQAGILPSGEPSSAEGSFAVAWPSQTAEPGPAQPAEASEVVGGAQEPGETAASEATSSSLPAVVVETFSATVNGAVEGSAGTGRLDLPPGFMFKVQAQHDYTATDTDELQLKAGDVVLVIPFQNPEEQDEGWLMGVKESDWNQHKELEKCRGVFPENFTERVQ"
peptide = string[299:308]
print peptide
accession = "O08539"
phosPos = 304
print extractPeptide(peptide, accession, string, phosPos, 9)
