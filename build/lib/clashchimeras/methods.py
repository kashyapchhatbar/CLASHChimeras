import re

from itertools import repeat


def convertCigar(cigarString):
    cigarAttributes = re.findall("\D", cigarString)
    cigarNumbers = re.findall("\d+", cigarString)
    convertedCigarString = ""
    for index, attribute in enumerate(cigarAttributes):
        for i in repeat(attribute, int(cigarNumbers[index])):
            convertedCigarString += i
    return convertedCigarString


def chimeraOrNot(bitOne, bitTwo, overlap=4, gap=9):
    combinedCigar = ""
    for i, j in zip(bitOne, bitTwo):
        if i == "M" and j == "M":
            combinedCigar += "="
        elif i == "S" and j == "S":
            combinedCigar += "#"
        elif i == "I" and j == "D":
            combinedCigar += "-"
        elif i == "D" and j == "I":
            combinedCigar += "-"
        elif i == "M" and j != "M":
            combinedCigar += "{"
        elif i != "M" and j == "M":
            combinedCigar += "}"
        elif i != "M" and j == "D":
            combinedCigar += "-"
        elif i == "D" and j != "M":
            combinedCigar += "-"
        elif i == "I" and j != "M":
            combinedCigar += "+"
        elif i != "M" and j == "I":
            combinedCigar += "+"
        else:
            combinedCigar += "*"
    numberOfMs = len(list(filter(lambda x: x == "=", combinedCigar)))
    curlyStartStart = combinedCigar.find("{")
    curlyStartEnd = combinedCigar.rfind("{")
    curlyEndStart = combinedCigar.find("}")
    curlyEndEnd = combinedCigar.rfind("}")
    listCurlies = sorted([[curlyStartStart, curlyStartEnd], [curlyEndStart,
                                                             curlyEndEnd]])
    matches = combinedCigar.count("{") + combinedCigar.count("}") + \
        combinedCigar.count("=") + combinedCigar.count("-") + \
        combinedCigar.count("+")
    if numberOfMs <= overlap and abs(listCurlies[0][1] - listCurlies[1][0]) \
            <= gap and matches >= (0.75 * len(combinedCigar)):
        return combinedCigar
    else:
        return None


def findRegion(record_data):
    if 'ENS' in record_data.refId and '|' in record_data.refId:
        reference_id = record_data.refId
        reference_start = int(record_data.start)
        match_length = record_data.matchLength
        regions = reference_id.split("|")[7:-1]
        match_start = reference_start
        match_end = reference_start + match_length
        mixed = []

        for region in regions:
            region_name = region.split(":")[0]
            region_start = int(region.split(":")[-1].split("-")[0])
            region_end = int(region.split(":")[-1].split("-")[-1])
            if match_start >= region_start and match_end <= region_end:
                mixed.append(region_name)
                break
            if match_start >= region_start and match_start <= region_end:
                mixed.append(region_name)
            if match_end >= region_start and match_end <= region_end:
                mixed.append(region_name)
            if match_start <= region_start and match_end >= region_end:
                mixed.append(region_name)
        return '{}|{}|{}'.format(reference_id.split("|")[0],
                                 reference_id.split("|")[5],
                                 "-".join(mixed))
    else:
        return record_data.refId
