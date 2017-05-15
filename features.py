from sequence import sequence
from CTD import CalculateCTD
from QuasiSequenceOrder import GetQuasiSequenceOrder, GetSequenceOrderCouplingNumberTotal
import math

AA_1_letter = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

''' Constant '''
accession_index = 0
sequence_index = 1
position_index = 2

''' Feature 1: Shannon Entropy '''


def entropy_Window(string):
    frequency = {key: 0 for key in AA_1_letter}
    for i in range(len(string)):
        if string[i] != 'X':
            frequency[string[i]] += 1

    probability = {key: (value / float(len(string))) for key, value in frequency.items()}
    entropyElements = []
    for item in probability:
        if probability[item] > 0:
            entropyElements.append(-probability[item] * math.log(probability[item], 2))

    return sum(entropyElements)


def feature_1_shannon_entropy(sequence, window=9):
    phosPos = sequence.getPhosphorylationPosition()

    sequenceString = sequence.getValue()
    n = sequence.getSize()

    # frequency = {key: 0 for key in AA_1_letter}
    # for i in range(n):
    #     frequency[sequenceString[i]] += 1
    #
    # probability = {key: (value / float(n)) for key, value in frequency.items()}
    # entropyElements = []
    # for item in probability:
    #     if probability[item] > 0:
    #         entropyElements.append(-probability[item] * math.log(probability[item], 2))

    if phosPos < window / 2:
        return entropy_Window(sequenceString[:phosPos - 1 + window / 2 + 1])

    if phosPos >= n - window / 2:
        return entropy_Window(sequenceString[phosPos - 1 - window / 2:])

    return entropy_Window(sequenceString[phosPos - window / 2:phosPos + window / 2 + 1])
    # return sum(entropyElements)

    # probability = {}
    # sequenceString = sequence.getValue()
    #
    # i = -1
    # # print sequenceString
    # while i < len(sequenceString):
    #     i += 1
    #     if i % 2 == 0:
    #         if sequenceString[i] in probability:
    #             probability[sequenceString[i]] += 1
    #         else:
    #             probability[sequenceString[i]] = 1
    #
    # probability = {probability[key]/float(len(sequenceString)/2+1) for key in probability}
    #
    # for key in probability:
    #     probability[key] /= float(len(sequenceString) / 2 + 1)
    #
    # entropy = [-probability[key] * math.log(probability[key], 2) for key in probability]
    # return sum(entropy)


''' Feature 2: Relative Entropy '''


def relativeEntropy_Window(string):
    frequency = {key: 0 for key in AA_1_letter}
    uniDis = 1 / float(len(string))

    for i in range(len(string)):
        if string[i] != 'X':
            frequency[string[i]] += 1

    probability = {key: (value / float(len(string))) for key, value in frequency.items()}
    KLDistance = []
    for item in probability:
        if probability[item] > 0:
            KLDistance.append(probability[item] * math.log(probability[item] / uniDis, 2))

    return sum(KLDistance)


def feature_2_relative_entropy(sequence, window=9):
    phosPos = sequence.getPhosphorylationPosition()

    sequenceString = sequence.getValue()
    n = sequence.getSize()

    # frequency = {key: 0 for key in AA_1_letter}
    # for i in range(n):
    #     frequency[sequenceString[i]] += 1
    #
    # probability = {key: (value / float(n)) for key, value in frequency.items()}
    #
    # KLdistance = {}
    # uniDis = 1/float(len(sequenceString))
    # for key in probability:
    #     if probability[key] != 0:
    #         KLdistance[key] = probability[key] * math.log(probability[key] / uniDis, 2)
    #     else:
    #         KLdistance[key] = 0

    if phosPos < window / 2:
        return relativeEntropy_Window(sequenceString[:phosPos + window / 2 + 1])

    if phosPos >= n - window / 2:
        return relativeEntropy_Window(sequenceString[phosPos - window / 2:])

    return relativeEntropy_Window(sequenceString[phosPos - window / 2:phosPos + window / 2 + 1])
    # return sum(KLdistance.values())


''' Feature 3: Information Gain'''


def feature_3_information_gain(sequence, window=9):
    """ Information Gain: IG = H - RE """
    return feature_1_shannon_entropy(sequence, window) - feature_2_relative_entropy(sequence, window)


''' Feature 4-12: Solvent Accessible Surface (ASA)'''


# TODO

def feature_4_12_ASA(sequence, window=9):
    sequenceString = sequence.getValue()
    n = sequence.getSize()

    position = sequence.getPhosphorylationPosition()

    accession = sequence.getAccession()
    phosPos = sequence.getPhosphorylationPosition()

    if phosPos < window / 2:
        windowString = sequenceString[:phosPos + window / 2 + 1]
        begin = 0
        end = phosPos + window / 2
    else:
        if phosPos >= n - window / 2:
            windowString = sequenceString[phosPos - window / 2:]
            begin = phosPos - window / 2
            end = n - 1
        else:
            windowString = sequenceString[phosPos - window / 2:phosPos + window / 2 + 1]
            begin = phosPos - window / 2
            end = phosPos + window / 2

    f_asa = open("./RVP_seq_output/" + accession + ".out")
    tmp_list = []
    while True:
        strTmp = f_asa.readline()
        if strTmp == "":
            break

        tmp_list.append(strTmp)

    tmp_list = tmp_list[int(position) - 1 - window / 2: int(position) - 1 + window / 2 + 1]
    tmp_list = [item.split("\t") for item in tmp_list]

    tmp_list = [item[1] for item in tmp_list]
    tmp_list = [item[1:len(item) - 3] for item in tmp_list]

    return tmp_list


''' Feature 13-102: Overlapping properties'''


def feature_13_102_overlapping_properties(sequence):
    sequenceString = sequence.getWindowValue()
    # print sequenceString

    polar = 'NQSDECTKRHYW'
    positive = 'KHR'
    negative = 'DE'
    charged = 'KHRDE'
    hydrophobic = 'AGCTIVLKHFWYM'
    aliphatic = 'IVL'
    aromatic = 'FYWH'
    small = 'PNDTCAGSV'
    tiny = 'ASGC'
    proline = 'P'

    def aa_overlapping_properties(aa):
        result = ''
        result += '1,' if aa in polar else '0,'
        result += '1,' if aa in positive else '0,'
        result += '1,' if aa in negative else '0,'
        result += '1,' if aa in charged else '0,'
        result += '1,' if aa in hydrophobic else '0,'
        result += '1,' if aa in aliphatic else '0,'
        result += '1,' if aa in aromatic else '0,'
        result += '1,' if aa in small else '0,'
        result += '1,' if aa in tiny else '0,'
        result += '1,' if aa in proline else '0,'

        return result[:len(result) - 1]

    sequence_result = ''
    for aa_index in range(len(sequenceString)):
        sequence_result += aa_overlapping_properties(sequenceString[aa_index])
        sequence_result += ','

    return sequence_result[:len(sequence_result) - 1]


''' Feature 103-106: Average Cumulative Hydrophobicity '''


def feature_103_106_hydrophobicity(sequence):
    hydrophobicity_scales = sequence.getHydrophobicity()

    def generateSubwindowssequence(sequence):
        # KAGVSPHED
        sequenceString = sequence.getWindowValue()
        size = len(sequenceString)

        center = size / 2 + 1  # center = 5
        numberOfFeatures = size - center  # num = 4

        tmp = []
        for i in range(1, numberOfFeatures + 1):
            tmp.append(sequenceString[center - i - 1:center + i])

        return tmp

    subWindows = generateSubwindowssequence(sequence)

    f = {}  # feature array
    for item in subWindows:
        sum = 0
        for i in range(len(item)):
            if item[i] != "X":
                sum += hydrophobicity_scales[item[i]]

        f[item] = sum / float(len(item))
    return f.values()


''' Feature 107-286: Sequence Feature '''
#  This part is intentionally skipped
#  [22], springer

''' Feature 287-433: Composition, Transition, and Distribution (CTD) '''


# Old code moved to archived.py
# This part used code from propy package

def feature_287_433_CTD(sequence):
    return CalculateCTD(sequence.getValue())


''' Feature 434-493: Sequence Order Coupling Numbers '''


# Old code moved to archived.py
# This part used code from propy package
#  PROBLEM: + maxlag
#           + lack of elements (the OP doesn't calculating on Grantham distance matrix)

def feature_434_493_socn(sequence, maxlag=30):
    phosPos = sequence.getPhosphorylationPosition()
    sequenceString = sequence.getValue()

    if maxlag % 2 == 0:
        window = maxlag + 1
    else:
        window = maxlag

    n = sequence.getSize()

    if phosPos < window / 2:
        return GetSequenceOrderCouplingNumberTotal(sequenceString[:phosPos + window / 2 + 1])

    if phosPos >= n - window / 2:
        return GetSequenceOrderCouplingNumberTotal(sequenceString[phosPos - window / 2:])

    return GetSequenceOrderCouplingNumberTotal(sequenceString[phosPos - window / 2:phosPos + window / 2 + 1])

    # return GetSequenceOrderCouplingNumberTotal(sequence.getValue())


''' Feature 494-593: Quasi Sequence Order (QSO)'''


# Old code moved to archived.py
# This part used code from propy package

def feature_494_593_quasi(sequence, maxlag=30):
    phosPos = sequence.getPhosphorylationPosition()
    sequenceString = sequence.getValue()

    if maxlag % 2 == 0:
        window = maxlag + 1
    else:
        window = maxlag

    n = sequence.getSize()

    if phosPos < window / 2:
        return GetQuasiSequenceOrder(sequenceString[:phosPos + window / 2 + 1])

    if phosPos >= n - window / 2:
        return GetQuasiSequenceOrder(sequenceString[phosPos - window / 2:])

    return GetQuasiSequenceOrder(sequenceString[phosPos - window / 2:phosPos + window / 2 + 1])
    # return GetQuasiSequenceOrder(sequence.getValue())

# Test part
# test = feature_1_shannon_entropy3_102(sequence('KAGVSPHED'))
# test = feature_1_shannon_entropy03_106(sequence('KAGVSPHED'))
# test = feature_2_relative_entropy87_433(sequence('RKEDQNGASTPHYCLVIMFW'))
# testString = sequence('ELRLRYCAPAGFALLKCNDADYDGFKTNCSNVSVVHCTNLMNTTVTTGLLLNGSYSENRTQIWQKHRTSNDSALILLNKHYNLTVTCKRPGNKTVLPVTIMAGLVFHSQKYNLRLRQAWCHFPSNWKGAWKEVKEEIVNLPKERYRGTNDPKRIFFQRQWGDPETANLWFNCHGEFFYCKMDWFLNYLNNLTVDADHNECKNTSGTKSGNKRAPGPCVQRTYVACHIRSVIIWLETISKKTYAPPREGHLECTSTVTGMTVELNYIPKNRTNVTLSPQIESIWAAELDRYKLVEITPIGFAPTEVRRYTGGHERQKRVPFVVQSQHLLAGILQQQKNLLAAVEAQQQMLKLTIWGVK')
# testString = sequence('ELDDQRTYR', 5)
#
# test1 = feature_1_shannon_entropy(testString)
# test2 = feature_2_relative_entropy(testString)
# test3 = feature_3_information_gain(testString)
# test4 = feature_13_102_overlapping_properties(testString)
# test5 = feature_103_106_hydrophobicity(testString)
# test6 = feature_287_433_CTD(testString)
# test7 = feature_434_493_socn(testString, 4)
# test8 = feature_494_593_quasi(testString)
#
# print test1
# print test2
# print test3
# print test4
# print test5
# print test6
# print test7
# print test8
