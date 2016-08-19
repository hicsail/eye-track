import numpy as np
import csv
from datetime import datetime
from datetime import timedelta
import math
import random
import time
import heapq

#xp = [0, 1,2,3]#x-values
#fp = [-3, 99.27,2,0]#y-values
#print(np.interp(2.5, xp, fp))#interpolated x-coordinates, x-values, y-values'''

#modInterp is basically just interp, except if you try to interpolate data outside of the given range of x-values,
#it still works
def modInterp(predX, xvals, yvals):
    if predX > max(xvals):
        p1 = (max(xvals),yvals[xvals.index(max(xvals))])
        secondBiggest = heapq.nlargest(2,xvals)[-1]
        p2 = (secondBiggest,yvals[xvals.index(secondBiggest)])
        m = float(p1[1]-p2[1])/float(p1[0]-p2[0])
        b = p1[1] - m*(p1[0])
        return (m*predX+b)
    elif predX < min(xvals):
        p1 = (min(xvals),yvals[xvals.index(min(xvals))])
        secondSmallest = heapq.nsmallest(2,xvals)[-1]
        p2 = (secondSmallest,yvals[xvals.index(secondSmallest)])
        m = float(p1[1]-p2[1])/float(p1[0]-p2[0])
        b = p1[1] - m*(p1[0])
        return (m*predX+b)
    else:
        return (np.interp(predX, xvals, yvals))
#print(modInterp(-1, [0,1,2,3], [-3,99.27,2,0]))
#print(modInterp(3.5, [0,1,2,3], [-3,99.27,2,0]))

#Option to have original non-conforming entries "forgotten" when writing to the new file, as well as
#   option of FPS of which ms intervals to interpolate,
#   option of extrapolating data (for a given trial) outside the range of entries based on the last 2 data entries

#She said that the lastData and delete options should be, by default, True
def interpolate(filename, outputFile, fps, delete=True, lastData=True):
    #begin = time.time()
    if (fps <= 0):
        print("Cannot have a non-positive value for fps")
        return

    #turning fps into milliseconds-per-frame
    mspf = (1.0/fps)*1000

    with open(filename) as tsv, open(outputFile, 'w') as out:
        read = csv.reader(tsv, delimiter = "\t")

        #write the first line of the old file into the new format, which
        #are the headers/categories
        firstLine = next(read)
        out.write("\t".join(firstLine)+"\n")

        #startMS is just to keep up with the updated RecordingTime [ms]
        #startLine is the first line of each "Trial", AKA the first actual entry
        #   after each 'Separator' line between trials
        #prev is the line before the current read line, for interpolation purposes (since
        #   we're only using the last data-entry and the current data-entry to interpolate
        #   a new entry)
        #prev2 is for the purposes of the option 'lastData', in that it's the line before
        #   the previous line (since at the end of a trial's data range, we're using the
        #   last 2 data-entries to interpolate, and you only know that you're at the end
        #   of the previous trial because you're at the beginning 'Separator' line of the next trial)
        startMS = 0
        startLine = []
        prev = []
        prev2 = []
        #This was just to parse through the first 39 or so lines of the file for testing purposes, rather than
        #   the entire file
        #i = 0
        for line in read:
            #boolean condition to come into place for the last 'else'-condition
            smallGap = False
            while not smallGap:
                if line[10] == 'Separator':#At the start of a new trial line, demarcated by the 'Separator' in the 'Category'
                    #If the option of generating the last-data entry for a given trial is marked true, before updating
                    #the startLine and prev-line references, write to the new file this new interpolated entry
                    if lastData:
                        if startMS > 0:#Just so that the first iteration is skipped, since startMS is initialized to 0 outside the loop
                            #Update RecordingTime [ms] with mspf since this interpolated entry's timestamp is just the last timestamp + mspf
                            startMS += mspf

                            #Convert the 'Time of Day' category of the previous line into a Python datetime object, then add to the date
                            #the difference in milliseconds between the interpolated entry's RecordingTime and the previous entry's RecordingTime
                            timeOfDay = datetime.strptime(prev[1] + "000", "%H:%M:%S:%f") + timedelta(milliseconds = (float(startMS) - float(prev[0])))

                            #The following commented-out code was just for rounding purposes -- but I assumed the default was to always round-down
                            #the remaining microseconds
                            #if (timeOfDay.microsecond/1000.0) - int(timeOfDay.microsecond/1000.0) < 0.5:
                            #    timeOfDay = timeOfDay.strftime('%H:%M:%S:%f')[:-3]
                            #else:
                            #    timeOfDay = (timeOfDay + timedelta(milliseconds=1)).strftime('%H:%M:%S:%f')[:-3]

                            #Turn timeOfDay into a string, ignoring the last 3 digits of the microsecond-places
                            timeOfDay = timeOfDay.strftime('%H:%M:%S:%f')[:-3]

                            #Basically just copy the previous line's filled categories, then modify the following index-categories
                            newLine = list(prev)
                            newLine[0] = str(round(startMS, 1))#round to one-place to the one-tenth-millisecond
                            newLine[1] = timeOfDay
                            #Interpolate the following categories, using last 2 data-entries, which would be stored in the prev and prev2 references
                            #Index 12: Pupil Diameter Right
                            newLine[12] = str(modInterp(startMS, [float(prev[0]), float(prev2[0])], [float(prev[12]), float(prev2[12])]))
                            #Index 13: Point of Regard Right X
                            newLine[13] = str(modInterp(startMS, [float(prev[0]), float(prev2[0])], [float(prev[13]), float(prev2[13])]))
                            #Index 14: Point of Regard Right Y
                            newLine[14] = str(modInterp(startMS, [float(prev[0]), float(prev2[0])], [float(prev[14]), float(prev2[14])]))
                            #Index 16: Gaze Vector Right X
                            newLine[16] = str(modInterp(startMS, [float(prev[0]), float(prev2[0])], [float(prev[16]), float(prev2[16])]))
                            #Index 17: Gaze Vector Right Y
                            newLine[17] = str(modInterp(startMS, [float(prev[0]), float(prev2[0])], [float(prev[17]), float(prev2[17])]))
                            #Index 18: Gaze Vector Right Z
                            newLine[18] = str(modInterp(startMS, [float(prev[0]), float(prev2[0])], [float(prev[18]), float(prev2[18])]))

                            #Comment out later -- was just to see which entries in the new file were interpolated
                            #newLine.append('Interpolated')
                            
                            out.write("\t".join(newLine)+"\n")

                    #Write the 'Separator' line, marking the start of the trial, to the new file
		    #-->She said she didn't want the 'Separator' line in the final output file
                    #out.write("\t".join(line)+"\n")
                    
                    #Initialize the prev and startLine references to the first entry after the 'Separator' line,
                    #as well as the startMS for the beginning RecordingTime of the trial
                    startLine = next(read)
                    prev = startLine
                    startMS = float(startLine[0])
                    out.write("\t".join(startLine)+"\n")
                    smallGap = True
                #If the current read-line's RecordingTime - the startMS of the current time interval < mspf:
                elif (float(line[0]) - startMS) < mspf:
                    #If the option of lastData is true, also update the prev2 reference
                    if lastData:
                        prev2 = prev
                    #Update prev reference to the current line, for the next iteration of the for-loop
                    prev = line
                    #If we ARE writing the old data-entries, non-conforming with the constant frame-rate, write them
                    if not delete:
                        out.write("\t".join(line)+"\n")
                    smallGap = True
                #If the current read-line's RecordingTime - the startMS of the current time interval == mspf,
                #AKA, the current data-entry conforms with the specified constant frame-rate -- write the current line
                elif (float(line[0]) - startMS) == mspf:
                    if lastData:
                        prev2 = prev
                    #The new time interval begins with the RecordingTime of this data-entry
                    startMS = float(line[0])
                    prev = line
                    startLine = line
                    out.write("\t".join(line)+"\n")
                    smallGap = True
                #Else, if the current read-line's RecordingTime - the startMS of the current time interval < mspf,
                #AKA, we need to interpolate a new data-entry, create and write the entry:
                else:
                    #Update the startMS to the RecordingTime of this interpolated data-entry, using the last startMS -- the beginning
                    #of the last time-interval -- and adding mspf
                    startMS += mspf

                    #Convert the 'Time of Day' category of the previous line into a Python datetime object, then add to the date
                    #the difference in milliseconds between the interpolated entry's RecordingTime and the previous entry's RecordingTime
                    timeOfDay = datetime.strptime(prev[1] + "000", "%H:%M:%S:%f") + timedelta(milliseconds = (float(startMS) - float(prev[0])))

                    #Option of rounding time up or down
                    #if (timeOfDay.microsecond/1000.0) - int(timeOfDay.microsecond/1000.0) < 0.5:
                    #    timeOfDay = timeOfDay.strftime('%H:%M:%S:%f')[:-3]
                    #else:
                    #    timeOfDay = (timeOfDay + timedelta(milliseconds=1)).strftime('%H:%M:%S:%f')[:-3]

                    #...But I went with the default of just always rounding down
                    #Turning timeOfDay into a string with the last 3 microsecond-places left out
                    timeOfDay = timeOfDay.strftime('%H:%M:%S:%f')[:-3]

                    #Basically, for the non-interpolated categories, we have the option of using the previous line's filled-in data
                    #or the current line's filled-in data -- so, on Frederick's suggestion, whichever data-entry is closer in time to
                    #the interpolated entry's RecordingTime gets precedence -- and if they're equally distant in time,
                    #idk just use Math.Random and flip a coin
                    newLine = []
                    if (startMS - float(prev[0])) < (float(line[0]) - startMS):
                        newLine = list(prev)
                    elif (startMS - float(prev[0])) > (float(line[0]) - startMS):
                        newLine = list(line)
                    else:#On Frederick's suggestion...
                        if random.random() > 0.5:
                            newLine = list(prev)
                        else:
                            newLine = list(line)

                    #Modify the following index-categories
                    newLine[0] = str(round(startMS, 1))
                    newLine[1] = timeOfDay
                    #Index 12: Pupil Diameter Right
                    newLine[12] = str(np.interp(startMS, [float(prev[0]), float(line[0])], [float(prev[12]), float(line[12])]))
                    #Index 13: Point of Regard Right X
                    newLine[13] = str(np.interp(startMS, [float(prev[0]), float(line[0])], [float(prev[13]), float(line[13])]))
                    #Index 14: Point of Regard Right Y
                    newLine[14] = str(np.interp(startMS, [float(prev[0]), float(line[0])], [float(prev[14]), float(line[14])]))
                    #Index 16: Gaze Vector Right X
                    newLine[16] = str(np.interp(startMS, [float(prev[0]), float(line[0])], [float(prev[16]), float(line[16])]))
                    #Index 17: Gaze Vector Right Y
                    newLine[17] = str(np.interp(startMS, [float(prev[0]), float(line[0])], [float(prev[17]), float(line[17])]))
                    #Index 18: Gaze Vector Right Z
                    newLine[18] = str(np.interp(startMS, [float(prev[0]), float(line[0])], [float(prev[18]), float(line[18])]))

                    #Comment out later -- was just to see which entries in the new file were interpolated
                    #newLine.append('Interpolated')
                    
                    out.write("\t".join(newLine)+"\n")

                    #Update startLine, prev2, and prev references for the next iterations of the for-loop
                    startLine = newLine
                    if lastData:
                        prev2 = prev
                    prev = newLine

                    #This is where the while-loop condition comes into place;
                    #if the difference between the new interpolated entry's RecordingTime and the current line's RecordingTime
                    #is STILL larger than the specified frame-rate interval, then we have to interpolate more data-entries
                    #to fill in this 'large gap' of milliseconds
                    if (float(line[0]) - startMS) < mspf:
                        smallGap = True
                        if lastData:
                            prev2 = prev
                        prev = line
                        if not delete:
                            out.write("\t".join(line)+"\n")
            #if i == 39:
            #    break
            #i += 1

        #Just to see how long the whole script takes
        #print("Finished in {0} seconds".format(time.time() - begin))

#Over 28,497 entries: takes ~2.7 seconds
interpolate('new.tsv', 'newOutput.tsv', 75, delete=True)




















