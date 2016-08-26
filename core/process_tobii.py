import csv

#Slight implementation note: I actually accidentally wrote all the string-values to the outputfile.csv as 'some_string', whereas
#I didn't notice after the fact that all such string values in the original input files are written as '"some_string"', and written to
#the outputfile.csv that way...
def process_tobii(tobiifile, aoifile, phasefile, condition, trialorder, outputfile):
    tobii = []
    aoi = []
    phase = []
    #Read the files in as 2D arrays
    with open(tobiifile) as t, open(aoifile) as a, open(phasefile) as phaseTiming:
        tobii = [line.split("\t")[:-1] for line in t]
        aoi = [line.strip().split("\t") for line in a]
        phase = [line.strip().split("\t") for line in phaseTiming]

    #I declared a ton of variables that would make it easier (in case certain column-indices switched around in different file formats)
    #to simply store the given index of a specific column at the beginning (from the column-headers), and then continually use that stored index throughout,
    #by way of just a key into a dict
    tobiiIndices = {}
    aoiIndices = {}
    phaseIndices = {}

    #Create the additional headers for the created-columns now for the tobii-file
    tobii[0] += ["X", "Condition", "TrialOrder", "Trial", "Phase", "FramesFromMovieOnset", "FramesFromPhaseOnset", "TimeFromMovieOnset", "TimeFromPhaseOnset", 
		"Subphase", "FramesFromSubphaseOnset", "TimeFromSubphaseOnset", "GazeXAvg", "GazeYAvg", "SceneName", "Position", "SceneType"]
    for i in range(len(tobii[0])):
        #rename some columns
        if tobii[0][i] == "GazePointLeftX..ADCSpx." or tobii[0][i] == "GazePointLeftX (ADCSpx)":
            tobii[0][i] = "GazePointLeftX"
        elif tobii[0][i] == "GazePointRightX..ADCSpx." or tobii[0][i] == "GazePointRightX (ADCSpx)":
            tobii[0][i] = "GazePointRightX"
        elif tobii[0][i] == "GazePointLeftY..ADCSpx." or tobii[0][i] == "GazePointLeftY (ADCSpx)":
            tobii[0][i] = "GazePointLeftY"
        elif tobii[0][i] == "GazePointRightY..ADCSpx." or tobii[0][i] == "GazePointRightY (ADCSpx)":
            tobii[0][i] = "GazePointRightY"
        elif tobii[0][i] == "X.StudioTestName":
            tobii[0][i] = "StudioTest"
        elif tobii[0][i] == "X.gender.Value" or tobii[0][i] == "[gender]Value":
            tobii[0][i] = "Gender"
        elif tobii[0][i] == "MediaName":
            tobii[0][i] = "MovieName"
        #store the index of that column
        tobiiIndices[tobii[0][i]] = i
    #refer to them later easily by memorable variable-names
    movie = tobiiIndices["MovieName"]
    tobiiParticipant = tobiiIndices["ParticipantName"]
    event = tobiiIndices["StudioEvent"]
    gRY = tobiiIndices["GazePointRightY"]
    gRX = tobiiIndices["GazePointRightX"]
    gLY = tobiiIndices["GazePointLeftY"]
    gLX = tobiiIndices["GazePointLeftX"]
    vL = tobiiIndices["ValidityLeft"]
    vR = tobiiIndices["ValidityRight"]
    cond = tobiiIndices["Condition"]
    torder = tobiiIndices["TrialOrder"]
    tobiiTrial = tobiiIndices["Trial"]
    tobiiPhase = tobiiIndices["Phase"]
    fMovie = tobiiIndices["FramesFromMovieOnset"]
    fPhase = tobiiIndices["FramesFromPhaseOnset"]
    tMovie = tobiiIndices["TimeFromMovieOnset"]
    tPhase = tobiiIndices["TimeFromPhaseOnset"]
    tobiiSubphase = tobiiIndices["Subphase"]
    fSubphase = tobiiIndices["FramesFromSubphaseOnset"]
    tSubphase = tobiiIndices["TimeFromSubphaseOnset"]
    gazeX = tobiiIndices["GazeXAvg"]
    gazeY = tobiiIndices["GazeYAvg"]
    tobiiScene = tobiiIndices["SceneName"]
    tobiiPos = tobiiIndices["Position"]
    tobiiType = tobiiIndices["SceneType"]

    #same for aoi-file
    for i in range(len(aoi[0])):
        aoiIndices[aoi[0][i]] = i
    aoiTrial = aoiIndices["Trial"]
    scene = aoiIndices["SceneName"]
    pos = aoiIndices["Position"]
    aoiType = aoiIndices["SceneType"]
    maxAOI = aoiIndices["Top Right (x,y)"]
    minAOI = aoiIndices["Bottom Left (x,y)"]

    #This was to store the entries of the aoi-file by using the filename as the key, and the row's entry
    #as the value -- such that I wouldn't have to loop over the entire aoi-file for every row in the tobii-file
    aoiTrials = {}
    for i in range(1, len(aoi)):#skip the first line of column-headers
        trialkey = aoi[i][aoiTrial]
        if trialkey in aoiTrials:
            aoiTrials[trialkey] += [aoi[i]]
        else:
            aoiTrials[trialkey] = [aoi[i]]

    for i in range(len(phase[0])):
        phaseIndices[phase[0][i]] = i
    filename = phaseIndices['Filename']
    trial = phaseIndices['Trial']
    phases = phaseIndices['Phase']
    subphase = phaseIndices['Subphase']
    start = phaseIndices['Starting_time']

    #Same concept here with the phasefile as with the aoi-file, except I also store the (Subphase, Starting_time) as the value,
    #with the key being the 'Trial' category
    phaseFilenames = {}
    subphaseTimes = {}
    for i in range(1, len(phase)):#skip the first line of column-headers
        filekey = phase[i][filename]
        sublist = phase[i]
        toAdd = [sublist[trial], sublist[phases]]
        if filekey in phaseFilenames:
            phaseFilenames[filekey] += [toAdd]
        else:
            phaseFilenames[filekey] = [toAdd]

        trialkey = phase[i][trial]
        if trialkey in subphaseTimes:
            subphaseTimes[trialkey] += [(phase[i][subphase], phase[i][start])]
        else:
            subphaseTimes[trialkey] = [(phase[i][subphase], phase[i][start])]

    hierarchy = {}
    index = 1#skip the first line of column-headers
    while (index != len(tobii)):
        #Essentially a while-loop that iterates over the tobii-file, but the incremented indexes are kept consistent
        #as I continually pop any rows that have the following condition:
        if tobii[index][event] == "MovieStart" or tobii[index][event] == "MovieEnd":
            tobii.pop(index)
        else:
            #toAppend is the enhanced row-entry with the extra columns, all with default '' values for now
            toAppend = tobii[index] + ['', condition, trialorder] + ['']*14

            #Simultaneously fill in the Trial and Phase categories for each row, while filling in any leading NA-gaps
            curMovie = tobii[index][movie]
            if curMovie in phaseFilenames:
                #0 for middle-index because I think, for a given filename, the trial and phase are consistent throughout the rows
                toAppend[tobiiTrial] = phaseFilenames[curMovie][0][0]
                toAppend[tobiiPhase] = phaseFilenames[curMovie][0][1]
            elif index == 1:
                pass
            else:
                toAppend[tobiiTrial] = tobii[index-1][tobiiTrial]

            #Here, what I'm doing is filtering out "blocks" of rows by ParticipantName --> Condition --> MovieName/Trial --> Phase --> Subphase,
            #by simply creating a hiearchy of nested dictionaries that each row-entry assigns itself to, incrementing the counter or becoming
            #the first (so the 0th frame) of that "block" of rows --> This is simpler than going through ParticipantName --> ... --> MovieName/Trial,
            #ParticipantName --> ... --> Subphase, and calculating separately the 'FramesFromSubphaseOnset' in several different for-loops
            curPartic = toAppend[tobiiParticipant]
            curTrial = toAppend[tobiiTrial]
            curPhase = toAppend[tobiiPhase]
            if curPartic != '' and condition != '' and curMovie != '':
                if curPartic not in hierarchy:
                    hierarchy[curPartic] = {condition: {curMovie: 0}}
                else:
                    if condition not in hierarchy[curPartic]:
                        hierarchy[curPartic][condition] = {curMovie: 0}
                    else:
                        if curMovie not in hierarchy[curPartic][condition]:
                            hierarchy[curPartic][condition][curMovie] = 0
                        else:
                            hierarchy[curPartic][condition][curMovie] += 1
                #Write to the 'FramesFromMovieOnset' and 'Time...' categories
                toAppend[fMovie] = hierarchy[curPartic][condition][curMovie]
                toAppend[tMovie] = round(toAppend[fMovie]*1000.0/60.0) #it just says "The data is collected at 60Hz"
                
                if curTrial != '' and curPhase != '':
                    if curTrial not in hierarchy[curPartic][condition]:
                        hierarchy[curPartic][condition][curTrial] = {curPhase: {'framesFromPhaseOnset': 0}}
                    else:
                        if curPhase not in hierarchy[curPartic][condition][curTrial]:
                            hierarchy[curPartic][condition][curTrial][curPhase] = {'framesFromPhaseOnset': 0}
                        else:
                            hierarchy[curPartic][condition][curTrial][curPhase]['framesFromPhaseOnset'] += 1
                    #Write to the 'FramesFromPhaseOnset' and 'Time...' categories
                    toAppend[fPhase] = hierarchy[curPartic][condition][curTrial][curPhase]['framesFromPhaseOnset']
                    toAppend[tPhase] = round(toAppend[fPhase]*1000.0/60.0) #it just says "The data is collected at 60Hz"

                    if curTrial in subphaseTimes:
                        #I iterate, for each row with a Subphase, through the a subset of entries in the phasefile to
                        #see the closest Starting_time it comes after
                        curTimestamp = toAppend[tPhase]
                        matchedSubphases = subphaseTimes[curTrial]
                        closestTime = 0
                        for x in range(1, len(matchedSubphases)):
                            closestDiff = curTimestamp - float(matchedSubphases[closestTime][1])*1000.0
                            curDiff = curTimestamp - float(matchedSubphases[x][1])*1000.0
                            if curDiff >= 0.0 and curDiff <= closestDiff:
                                closestTime = x
                        toAppend[tobiiSubphase] = matchedSubphases[closestTime][0]

                        curSubphase = toAppend[tobiiSubphase]
                        if curSubphase not in hierarchy[curPartic][condition][curTrial][curPhase]:
                            hierarchy[curPartic][condition][curTrial][curPhase][curSubphase] = 0
                        else:
                            hierarchy[curPartic][condition][curTrial][curPhase][curSubphase] += 1
                        #Write to the 'FramesFromSubphaseOnset' and 'Time...' categories
                        toAppend[fSubphase] = hierarchy[curPartic][condition][curTrial][curPhase][curSubphase]
                        toAppend[tSubphase] = round(toAppend[fSubphase]*1000.0/60.0) #it just says "The data is collected at 60Hz"
                else:
                    #These are just to stay consistent with how the original R-script wrote 0 as the default instead of NA, no matter what
                    toAppend[fSubphase] = 0
                    toAppend[tSubphase] = 0
            else:
                toAppend[fMovie] = 0
                toAppend[tMovie] = 0
                toAppend[fPhase] = 0
                toAppend[tPhase] = 0
                toAppend[fSubphase] = 0
                toAppend[tSubphase] = 0

            #Here, we compute the GazeX/YAvg
            cVL = float(toAppend[vL])
            cVR = float(toAppend[vR])
            #Some rows only had the right-eye GazePoint coordinates, some only had the left, some had neither -- so we only convert to float if they exist
            if toAppend[gRX] != '' and toAppend[gRY] != '':
                cRX = float(toAppend[gRX])
                cRY = float(toAppend[gRY])
            if toAppend[gLX] != '' and toAppend[gLY] != '':
                cLX = float(toAppend[gLX])
                cLY = float(toAppend[gLY])
            #Check validity and assign/compute GazeX/YAvg
            if (cVL > 1 or cVL < 0) and (cVR > 1 or cVR < 0):
                toAppend[gazeX] = -999
                toAppend[gazeY] = -999
            elif (cVL > 1 or cVL < 0) and (cVR <= 1 and cVR >= 0):
                #take the right side
                toAppend[gazeX] = cRX
                toAppend[gazeY] = cRY
            elif (cVR > 1 or cVR < 0) and (cVL <= 1 and cVL >= 0):
                #take the left side
                toAppend[gazeX] = cLX
                toAppend[gazeY] = cLY
            else:
                #avg both sides
                toAppend[gazeX] = (cLX + cRX)/2.0
                toAppend[gazeY] = (cLY + cRY)/2.0
            #We only attempt to write the SceneName, Position, and SceneType the following under these conditions:
            if toAppend[gazeX] != -999 and curTrial != '' and curTrial in aoiTrials:
                matchedAOI = aoiTrials[curTrial]
                for area in matchedAOI:
                    #Quite annoying -- parsing (x,y) coordinates where the numbers themselves can be written with commas if they're > 1000 -- like 45 lines of code for just this
                    maxCoords = area[maxAOI].strip()
                    minCoords = area[minAOI].strip()
                    maxCommas = maxCoords.count(',')
                    if maxCommas == 1:
                        maxXY = maxCoords.split(',')
                        xMax = int(maxXY[0])
                        yMax = int(maxXY[1])
                    elif maxCommas == 2:
                        maxXY = maxCoords.split(',')
                        temp1 = int(maxXY[0]) * (10 ** len(maxXY[1])) + int(maxXY[1])
                        temp2 = int(maxXY[1]) * (10 ** len(maxXY[2])) + int(maxXY[2])
                        if (temp1 - int(maxXY[2])) > (temp2 - int(maxXY[0])):
                            #use temp2
                            xMax = int(maxXY[0])
                            yMax = temp2
                        else:
                            #use temp1
                            xMax = temp1
                            yMax = int(maxXY[2])
                    else:#(For the case of 3 commas)
                        maxXY = maxCoords.split(',')
                        xMax = int(maxXY[0]) * (10 ** len(maxXY[1])) + int(maxXY[1])
                        yMax = int(maxXY[2]) * (10 ** len(maxXY[3])) + int(maxXY[3])
                    
                    minCommas = minCoords.count(',')
                    if minCommas == 1:
                        minXY = minCoords.split(',')
                        xMin = int(minXY[0])
                        yMin = int(minXY[1])
                    elif minCommas == 2:
                        minXY = minCoords.split(',')
                        temp1 = int(minXY[0]) * (10 ** len(minXY[1])) + int(minXY[1])
                        temp2 = int(minXY[1]) * (10 ** len(minXY[2])) + int(minXY[2])
                        if (temp1 - int(minXY[2])) > (temp2 - int(minXY[0])):
                            #use temp2
                            xMin = int(minXY[0])
                            yMin = temp2
                        else:
                            #use temp1
                            xMin = temp1
                            yMin = int(minXY[2])
                    else:#(For the case of 3 commas)
                        minXY = minCoords.split(',')
                        xMin = int(minXY[0]) * (10 ** len(minXY[1])) + int(minXY[1])
                        yMin = int(minXY[2]) * (10 ** len(minXY[3])) + int(minXY[3])
                    #Since I don't know the format around the pixel-coordinates of the corners of the area-of-interest on the eye-tracker screen, apparently
                    #the "Bottom-Left" y-pixel-coordinate value can be greater than the "Top-Right" y-pixel-coordinate value
                    if toAppend[gazeX] <= max(xMax, xMin) and toAppend[gazeX] >= min(xMax, xMin) and toAppend[gazeY] <= max(yMax, yMin) and toAppend[gazeY] >= min(yMax, yMin):
                        toAppend[tobiiScene] = area[scene]
                        toAppend[tobiiPos] = area[pos]
                        toAppend[tobiiType] = area[aoiType]
            #Otherwise, only in the case where ValidityLeft and ValidityRight were unacceptable, do we assign the value 'TrackLost''
            elif toAppend[gazeX] == -999:
                toAppend[tobiiScene] = "TrackLoss"
                toAppend[tobiiPos] = "TrackLoss"
                toAppend[tobiiType] = "TrackLoss"

            tobii[index] = toAppend
            index += 1
    #For some reason, only the 'StudioEvent' column could be left with '' values; all other '' turned into NA
    for i in range(len(tobii)):
        for j in range(len(tobii[i])):
            if j == event:
                continue
            elif tobii[i][j] == '':
                tobii[i][j] = 'NA'
            elif type(tobii[i][j]) != str:
                tobii[i][j] = str(tobii[i][j])#For numbers --> strings
    with open(outputfile, 'w') as out:
        for line in tobii:
            out.write(",".join(line)+"\n")#Comma-separated-values

process_tobii("oldinput_v2.tsv", "aoifile.tsv", "phasefile.tsv", "Intransitive", "Backward", "oldoutput_v2_python.csv")

#Checks whether a string represents a number
def RepresentsNumber(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False

def checkPython_process_tobii(R_outputfile, python_outputfile):
    R = []
    Py = []
    with open(R_outputfile) as Rout, open(python_outputfile) as Pout:
        R = [line.strip().split(",") for line in Rout]
        Py = [line.strip().split(",") for line in Pout]
    if (len(Py) != len(R)):
        print("Lengths aren't equal")
    else:
        for i in range(len(R)):
            for j in range(len(R[i])):
                if j < 26:
                    #columnChar is supposed to represent the first 52 column-headers Microsoft Excel displays for a .csv file (so A, B, C, ... AA, AB, AC, ...)
                    columnChar = str(chr(j+97)).upper()
                else:
                    columnChar = "A" + str(chr((j%26)+97)).upper()
                if i == 0 and j == 3:#This was the cell for 'Gender' -- somewhat strange how the original R-script changed the column-header
                    continue
                #If the value in the column is a string and not NA
                elif not RepresentsNumber(Py[i][j]) and Py[i][j] != "NA":
                    #Put back in the '"some_string"' double-quotes
                    if R[i][j] != ('\"' + Py[i][j] + '\"'):
                        print("Not equal at cell (1-indexed): ({0}, {1})".format(i+1, columnChar))
                        print("R-Output: {0}".format(R[i][j]))
                        print("Python-Output: {0}".format(Py[i][j]))
                        return
                elif RepresentsNumber(Py[i][j]):
                    #Cast both values in the R-output and the Python-output to float values, for precision purposes (as '4.10' != '4.1')
                    if float(R[i][j]) != float(Py[i][j]):
                        print("Not equal at cell (1-indexed): ({0}, {1})".format(i+1, columnChar))
                        print("R-Output: {0}".format(R[i][j]))
                        print("Python-Output: {0}".format(Py[i][j]))
                        return
                elif R[i][j] != Py[i][j]:
                    print("Not equal at cell (1-indexed): ({0}, {1})".format(i+1, columnChar))
                    print("R-Output: {0}".format(R[i][j]))
                    print("Python-Output: {0}".format(Py[i][j]))
                    return

checkPython_process_tobii("oldv2_output.csv", "oldoutput_v2_python.csv")



















