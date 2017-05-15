# # PROBLEM:
# # Composition: Group 1 and group 3 is replace by each other in the example (in the two papers)
# # RKEDQNGASTPHYCLVIMFW should be encoded as "33333322222221111111", but the paper show "11111122222223333333" encoding
# # Transition: Tmn should be 1/19 with 11111122222223333333 sequence (T12 = 1, T21 = 0)
# # This code will perform the feature as its definition, not the way the example shows us
# #  TODO
# def feature_287_433_CTD(sequence):
#     # Calculating 'Composition' first
#     encodedSequence = {}  # encoded sequence, with respect to attributes. It's calculated after calling compositionCalculation
#     n = sequence.getSize()  # Peptide size
#
#     CTD_group_7_features = sequence.getCTDGroup()  # Amino acids 7 attribute for grouping
#     CTD_attributes = sequence.getCTDattributes()  # CTD attributes list
#
#     sequenceString = sequence.getValue()  # Peptide
#
#     def compositionCalculation(sequence):
#
#         composition = {key: 0 for key in CTD_attributes}
#         for attribute in CTD_attributes:
#             result = ''
#             aa_per_group = {1: 0, 2: 0, 3: 0}
#
#             for i in range(len(sequenceString)):
#                 for groupIndex, CTD_group in enumerate(CTD_group_7_features):
#                     if sequenceString[i] in CTD_group[attribute]:
#                         result += str(groupIndex + 1)
#                         aa_per_group[groupIndex + 1] += 1
#                         break
#
#             encodedSequence[attribute] = result
#
#             composition[attribute] = {key: value / float(n) for key, value in aa_per_group.items()}
#
#         # return {key: value/float(n) for key, value in aa_per_group.items()}
#
#         # print encodedSequence
#         return composition
#
#     def transitionCalculation(sequence):
#         transition = {'12': 0, '13': 0, '23': 0}  # result of a transition of a seq with respect to a single attribute
#
#         transitionResult = {key: dict(transition) for key in CTD_attributes}
#
#         mn = ['12', '13', '23']  # class m to class n
#
#         for attribute in CTD_attributes:
#             for mnValue in mn:
#                 transitionResult[attribute][mnValue] = encodedSequence[attribute].count(mnValue) + encodedSequence[
#                     attribute].count(mnValue[::-1])
#         return transitionResult
#
#     compo = compositionCalculation(sequence)
#     trans = transitionCalculation(sequence)
#
#     return [compo, trans]











# def feature_434_493_socnmatrixName(sequenceVal, distanceMatrix, matrixName, maxlag):
#     tau = {}
#     n = sequenceVal.getSize()
#     sequenceString = sequenceVal.getValue()
#
#     for lag in range(1, maxlag + 1):
#         tmp_tau = 0
#         for i in range(n - lag):
#             # print i
#             tmp_tau += pow(distanceMatrix[sequenceString[i] + sequenceString[i + lag]], 2)
#
#         # tau_result['tau' + sequenceString[i] + sequenceString[i+lag]] = tmp_tau
#         # tau['tau ' + matrixName + str(lag)] = [lag, tmp_tau]
#         tau[lag] = ['tau ' + matrixName + str(lag), tmp_tau]
#     return tau
#
#
# # def feature_434_493_socnGran(sequenceVal, distanceMatrix, maxlag=30):
# #     tau = {}
# #     n = sequenceVal.getSize()
# #     sequenceString = sequenceVal.getValue()
# #
# #     for lag in range(1, maxlag + 1):
# #         tmp_tau = 0
# #         for i in range(n - lag):
# #             # print i
# #             tmp_tau += pow(distanceMatrix[sequenceString[i] + sequenceString[i + lag]], 2)
# #
# #         # tau_result['tau' + sequenceString[i] + sequenceString[i+lag]] = tmp_tau
# #         tau['tau Gran ' + str(lag)] = round(tmp_tau, 3)
# #         # print "n-lag" + str(n-lag)
# #     return tau
#
# def feature_434_493_socn(sequence, matrixName, maxlag=30):
#     # matrixName is either SW ( Schneider - Wrede ) or Gran (Grantham)
#     SW_matrix = sequence.getSchneiderWredeDistanceMatrix()
#     Gran_matrix = sequence.getGranthamDistanceMatrix()
#
#     tau_result = {}
#
#     if matrixName == 'SW':
#         tau_result.update(feature_434_493_socnmatrixName(sequence, SW_matrix, matrixName, maxlag))
#     if matrixName == 'Gran':
#         tau_result.update(feature_434_493_socnmatrixName(sequence, Gran_matrix, matrixName, maxlag))
#     return tau_result








# #  Literally probability of a letter to appear in the sequence
# def normalized_frequency(sequence):
#     sequenceString = sequence.getValue()
#     n = sequence.getSize()
#
#     frequency = {key: 0 for key in AA_1_letter}
#     for i in range(n):
#         frequency[sequenceString[i]] += 1
#
#     return {key: (value / float(n)) for key, value in frequency.items()}
#
#
# #  First 20 + 20 elements
# def feature_494_593_quasi_part1(sequence, weight, maxlag):
#     f = normalized_frequency(sequence)
#     f_sum = sum(f.values())  # f_sum = 1
#
#     quasi_1 = {}
#
#     socnSW = feature_434_493_socn(sequence, 'SW', maxlag)
#     socnGran = feature_434_493_socn(sequence, 'Gran', maxlag)
#
#     # tau_sum_SW = sum(feature_434_493(sequence, 'SW').values())
#     # tau_sum_Gran = sum(feature_434_493(sequence, 'Gran').values())
#
#     tau_sum_SW = sum([item[1] for item in socnSW.values()])
#     tau_sum_Gran = sum([item[1] for item in socnGran.values()])
#
#     denominatorSW = weight * tau_sum_SW + f_sum
#     denominatorGran = weight * tau_sum_Gran + f_sum
#
#     for aa in AA_1_letter:
#         numerator = f[aa]
#
#         #   Quasi
#         quasi_1['QSOSW ' + aa] = (numerator / float(denominatorSW))
#         quasi_1['QSOGran ' + aa] = (numerator / float(denominatorGran))
#
#     return quasi_1
#
#
# #  Last 30 + 30 elements
# def feature_494_593_quasi_part2(sequence, weight, maxlag=30):
#     f = normalized_frequency(sequence)
#     f_sum = sum(f.values())  # = 1
#
#     socnSW = feature_434_493_socn(sequence, 'SW', maxlag)
#     socnGran = feature_434_493_socn(sequence, 'Gran', maxlag)
#     denominatorSW = f_sum + weight * sum([item[1] for item in socnSW.values()])
#     denominatorGran = f_sum + weight * sum([item[1] for item in socnGran.values()])
#
#     result = {}
#     for index in range(21, 21 + maxlag):
#         # print "cac " + str(index) + " " + str(socnSW[index-20])
#         result['QSOSW' + str(index)] = round(weight * socnSW[index - 20][1] / denominatorSW, 4)
#         result['QSOGran' + str(index)] = round(weight * socnGran[index - 20][1] / denominatorGran, 4)
#
#     return result
#
#
# def feature_494_593_quasi(sequence, weight=0.1):
#     #   Others (80 elements)
#     quasiResult = {}
#     quasiResult.update(feature_494_593_quasi_part1(sequence, weight, maxlag=30))
#     quasiResult.update(feature_494_593_quasi_part2(sequence, weight, maxlag=30))
#
#     return quasiResult


def getASAList(file, position, window=9):
    tmp_list = []
    while True:
        strTmp = file.readline()
        if strTmp == "":
            break

        tmp_list.append(strTmp)

    tmp_list = tmp_list[int(position) - 1 - window / 2 : int(position) - 1 + window / 2 + 1]
    tmp_list = [item.split("\t") for item in tmp_list]

    tmp_list = [item[1] for item in tmp_list]
    tmp_list = [item[1:len(item)-3] for item in tmp_list]

    return tmp_list