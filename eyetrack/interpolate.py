import numpy as np
import csv
from datetime import datetime, timedelta
import random
# import time
import heapq
import argparse
import ast
from eyetrack import config


# mod_interpolate is basically just np.interp, except if you try to interpolate data outside of the given range of x-values,
# it still works
def mod_interpolate(predX, xvals, yvals):
    if predX > max(xvals):
        p1 = (max(xvals), yvals[xvals.index(max(xvals))])
        second_biggest = heapq.nlargest(2, xvals)[-1]
        p2 = (second_biggest, yvals[xvals.index(second_biggest)])
        m = float(p1[1] - p2[1]) / float(p1[0] - p2[0])
        b = p1[1] - m * (p1[0])
        return m * predX + b
    elif predX < min(xvals):
        p1 = (min(xvals), yvals[xvals.index(min(xvals))])
        second_smallest = heapq.nsmallest(2, xvals)[-1]
        p2 = (second_smallest, yvals[xvals.index(second_smallest)])
        m = float(p1[1] - p2[1]) / float(p1[0] - p2[0])
        b = p1[1] - m * (p1[0])
        return m * predX + b
    else:
        return np.interp(predX, xvals, yvals)


# Option to have original non-conforming entries "forgotten" when writing to the new file, as well as
#   option of FPS of which ms intervals to interpolate,
#   option of extrapolating data (for a given trial) outside the range of entries based on the last 2 data entries


def interpolate(input_file, output_file, fps, delete=True, last_data=True):
    # begin = time.time()
    if fps <= 0:
        print("Cannot have a non-positive value for fps")
        return

    # Turn fps into milliseconds-per-frame
    mspf = (1.0 / fps) * 1000

    with open(input_file, 'r', encoding='UTF-8') as tsv, open(output_file, 'w', encoding='UTF-8') as out:
        read = csv.reader(tsv, delimiter="\t")

        # write the first line of the old file into the new format, which
        # are the headers/categories
        first_line = next(read)

        # Column indexes and naming might change. By moving this to a config file, the user can set it themselves
        try:
            column_indexes = {'RecordingTime': first_line.index(config.get('RecordingTime')),
                              'TimeOfDay': first_line.index(config.get('TimeOfDay')),
                              'PupilRight': first_line.index(config.get('PupilRight')),
                              'RegardRightX': first_line.index(config.get('RegardRightX')),
                              'RegardRightY': first_line.index(config.get('RegardRightY')),
                              'GazeRightX': first_line.index(config.get('GazeRightX')),
                              'GazeRightY': first_line.index(config.get('GazeRightY')),
                              'GazeRightZ': first_line.index(config.get('GazeRightZ')),
                              'Category': first_line.index(config.get('Category'))}
        except ValueError as e:
            quit('Column name not found, ' + e.__str__())

        out.write("\t".join(first_line) + "\n")

        # startMS is just to keep up with the updated RecordingTime [ms]
        # startLine is the first line of each "Trial", AKA the first actual entry
        #   after each 'Separator' line between trials
        # prev is the line before the current read line, for interpolation purposes (since
        #   we're only using the last data-entry and the current data-entry to interpolate
        #   a new entry)
        # prev2 is for the purposes of the option 'lastData', in that it's the line before
        #   the previous line (since at the end of a trial's data range, we're using the
        #   last 2 data-entries to interpolate, and you only know that you're at the end
        #   of the previous trial because you're at the beginning 'Separator' line of the next trial)
        start_ms = 0
        start_line = []
        prev = []
        prev2 = []
        # This was just to parse through the first 39 or so lines of the file for testing purposes, rather than
        #   the entire file
        # i = 0
        for line in read:
            # boolean condition to come into place for the last 'else'-condition
            small_gap = False
            while not small_gap:
                # At the start of a new trial line, demarcated by the 'Separator' in the 'Category'
                if line[column_indexes['Category']] == 'Separator':
                    # If the option of generating the last-data entry for a given trial is marked true, before updating
                    # the startLine and prev-line references, write to the new file this new interpolated entry
                    if last_data:
                        # Just so that the first iteration is skipped, since startMS is initialized to 0 outside the loop
                        if start_ms > 0:
                            # Update RecordingTime [ms] with mspf since this interpolated entry's timestamp is just the last timestamp + mspf
                            start_ms += mspf

                            # Convert the 'Time of Day' category of the previous line into a Python datetime object, then add to the date
                            # the difference in milliseconds between the interpolated entry's RecordingTime and the previous entry's RecordingTime
                            time_of_day = datetime.strptime(prev[column_indexes['TimeOfDay']] + "000",
                                                            "%H:%M:%S:%f") + timedelta(
                                milliseconds=(float(start_ms) - float(prev[column_indexes['RecordingTime']])))

                            # The following commented-out code was just for rounding purposes -- but I assumed the default was to always round-down
                            # the remaining microseconds
                            # if (timeOfDay.microsecond/1000.0) - int(timeOfDay.microsecond/1000.0) < 0.5:
                            #    timeOfDay = timeOfDay.strftime('%H:%M:%S:%f')[:-3]
                            # else:
                            #    timeOfDay = (timeOfDay + timedelta(milliseconds=1)).strftime('%H:%M:%S:%f')[:-3]

                            # Turn timeOfDay into a string, ignoring the last 3 digits of the microsecond-places
                            time_of_day = time_of_day.strftime('%H:%M:%S:%f')[:-3]

                            # Basically just copy the previous line's filled categories, then modify the following index-categories
                            new_line = list(prev)
                            new_line[column_indexes['RecordingTime']] = str(
                                round(start_ms, 1))  # round to one-place to the one-tenth-millisecond
                            new_line[column_indexes['TimeOfDay']] = time_of_day
                            # Interpolate the following categories, using last 2 data-entries, which would be stored in the prev and prev2 references
                            # Index 12: Pupil Diameter Right
                            new_line[column_indexes['PupilRight']] = str(mod_interpolate(start_ms,
                                                                                         [float(prev[column_indexes[
                                                                                             'RecordingTime']]),
                                                                                          float(prev2[column_indexes[
                                                                                              'RecordingTime']])],
                                                                                         [float(prev[column_indexes[
                                                                                             'PupilRight']]),
                                                                                          float(prev2[column_indexes[
                                                                                              'PupilRight']])]))
                            # Index 13: Point of Regard Right X
                            new_line[column_indexes['RegardRightX']] = str(mod_interpolate(start_ms,
                                                                                           [float(prev[column_indexes[
                                                                                               'RecordingTime']]),
                                                                                            float(prev2[column_indexes[
                                                                                                'RecordingTime']])],
                                                                                           [float(prev[column_indexes[
                                                                                               'RegardRightX']]),
                                                                                            float(prev2[column_indexes[
                                                                                                'RegardRightX']])]))
                            # Index 14: Point of Regard Right Y
                            new_line[column_indexes['RegardRightY']] = str(mod_interpolate(start_ms,
                                                                                           [float(prev[column_indexes[
                                                                                               'RecordingTime']]),
                                                                                            float(prev2[column_indexes[
                                                                                                'RecordingTime']])],
                                                                                           [float(prev[column_indexes[
                                                                                               'RegardRightY']]),
                                                                                            float(prev2[column_indexes[
                                                                                                'RegardRightY']])]))
                            # Index 16: Gaze Vector Right X
                            new_line[column_indexes['GazeRightX']] = str(mod_interpolate(start_ms,
                                                                                         [float(prev[column_indexes[
                                                                                             'RecordingTime']]),
                                                                                          float(prev2[column_indexes[
                                                                                              'RecordingTime']])],
                                                                                         [float(prev[column_indexes[
                                                                                             'GazeRightX']]),
                                                                                          float(prev2[column_indexes[
                                                                                              'GazeRightX']])]))
                            # Index 17: Gaze Vector Right Y
                            new_line[column_indexes['GazeRightY']] = str(mod_interpolate(start_ms,
                                                                                         [float(prev[column_indexes[
                                                                                             'RecordingTime']]),
                                                                                          float(prev2[column_indexes[
                                                                                              'RecordingTime']])],
                                                                                         [float(prev[column_indexes[
                                                                                             'GazeRightY']]),
                                                                                          float(prev2[column_indexes[
                                                                                              'GazeRightY']])]))
                            # Index 18: Gaze Vector Right Z
                            new_line[column_indexes['GazeRightZ']] = str(mod_interpolate(start_ms,
                                                                                         [float(prev[column_indexes[
                                                                                             'RecordingTime']]),
                                                                                          float(prev2[column_indexes[
                                                                                              'RecordingTime']])],
                                                                                         [float(prev[column_indexes[
                                                                                             'GazeRightZ']]),
                                                                                          float(prev2[column_indexes[
                                                                                              'GazeRightZ']])]))

                            # Comment out later -- was just to see which entries in the new file were interpolated
                            # newLine.append('Interpolated')

                            out.write("\t".join(new_line) + "\n")

                    # Write the 'Separator' line, marking the start of the trial, to the new file
                    # -->She said she didn't want the 'Separator' line in the final output file
                    # out.write("\t".join(line)+"\n")

                    # Initialize the prev and startLine references to the first entry after the 'Separator' line,
                    # as well as the startMS for the beginning RecordingTime of the trial
                    start_line = next(read)
                    prev = start_line
                    start_ms = float(start_line[column_indexes['RecordingTime']])
                    out.write("\t".join(start_line) + "\n")
                    small_gap = True
                # If the current read-line's RecordingTime - the startMS of the current time interval < mspf:
                elif (float(line[column_indexes['RecordingTime']]) - start_ms) < mspf:
                    # If the option of lastData is true, also update the prev2 reference
                    if last_data:
                        prev2 = prev
                    # Update prev reference to the current line, for the next iteration of the for-loop
                    prev = line
                    # If we ARE writing the old data-entries, non-conforming with the constant frame-rate, write them
                    if not delete:
                        out.write("\t".join(line) + "\n")
                    small_gap = True
                # If the current read-line's RecordingTime - the startMS of the current time interval == mspf,
                # AKA, the current data-entry conforms with the specified constant frame-rate -- write the current line
                elif (float(line[column_indexes['RecordingTime']]) - start_ms) == mspf:
                    if last_data:
                        prev2 = prev
                    # The new time interval begins with the RecordingTime of this data-entry
                    start_ms = float(line[column_indexes['RecordingTime']])
                    prev = line
                    start_line = line
                    out.write("\t".join(line) + "\n")
                    small_gap = True
                # Else, if the current read-line's RecordingTime - the startMS of the current time interval < mspf,
                # AKA, we need to interpolate a new data-entry, create and write the entry:
                else:
                    # Update the startMS to the RecordingTime of this interpolated data-entry, using the last startMS -- the beginning
                    # of the last time-interval -- and adding mspf
                    start_ms += mspf

                    # Convert the 'Time of Day' category of the previous line into a Python datetime object, then add to the date
                    # the difference in milliseconds between the interpolated entry's RecordingTime and the previous entry's RecordingTime
                    time_of_day = datetime.strptime(prev[column_indexes['TimeOfDay']] + "000",
                                                    "%H:%M:%S:%f") + timedelta(
                        milliseconds=(float(start_ms) - float(prev[column_indexes['RecordingTime']])))

                    # Option of rounding time up or down
                    # if (timeOfDay.microsecond/1000.0) - int(timeOfDay.microsecond/1000.0) < 0.5:
                    #    timeOfDay = timeOfDay.strftime('%H:%M:%S:%f')[:-3]
                    # else:
                    #    timeOfDay = (timeOfDay + timedelta(milliseconds=1)).strftime('%H:%M:%S:%f')[:-3]

                    # ...But I went with the default of just always rounding down
                    # Turning timeOfDay into a string with the last 3 microsecond-places left out
                    time_of_day = time_of_day.strftime('%H:%M:%S:%f')[:-3]

                    # Basically, for the non-interpolated categories, we have the option of using the previous line's filled-in data
                    # or the current line's filled-in data -- so, on Frederick's suggestion, whichever data-entry is closer in time to
                    # the interpolated entry's RecordingTime gets precedence -- and if they're equally distant in time,
                    # idk just use Math.Random and flip a coin
                    if (start_ms - float(prev[column_indexes['RecordingTime']])) < (
                                float(line[column_indexes['RecordingTime']]) - start_ms):
                        new_line = list(prev)
                    elif (start_ms - float(prev[column_indexes['RecordingTime']])) > (
                                float(line[column_indexes['RecordingTime']]) - start_ms):
                        new_line = list(line)
                    else:  # On Frederick's suggestion...
                        if random.random() > 0.5:
                            new_line = list(prev)
                        else:
                            new_line = list(line)

                    # Modify the following index-categories
                    new_line[column_indexes['RecordingTime']] = str(round(start_ms, 1))
                    new_line[column_indexes['TimeOfDay']] = time_of_day
                    # Index 12: Pupil Diameter Right
                    new_line[column_indexes['PupilRight']] = str(
                        np.interp(start_ms, [float(prev[column_indexes['RecordingTime']]),
                                             float(line[column_indexes['RecordingTime']])],
                                  [float(prev[column_indexes['PupilRight']]),
                                   float(line[column_indexes['PupilRight']])]))
                    # Index 13: Point of Regard Right X
                    new_line[column_indexes['RegardRightX']] = str(
                        np.interp(start_ms, [float(prev[column_indexes['RecordingTime']]),
                                             float(line[column_indexes['RecordingTime']])],
                                  [float(prev[column_indexes['RegardRightX']]),
                                   float(line[column_indexes['RegardRightX']])]))
                    # Index 14: Point of Regard Right Y
                    new_line[column_indexes['RegardRightY']] = str(
                        np.interp(start_ms, [float(prev[column_indexes['RecordingTime']]),
                                             float(line[column_indexes['RecordingTime']])],
                                  [float(prev[column_indexes['RegardRightY']]),
                                   float(line[column_indexes['RegardRightY']])]))
                    # Index 16: Gaze Vector Right X
                    new_line[column_indexes['GazeRightX']] = str(
                        np.interp(start_ms, [float(prev[column_indexes['RecordingTime']]),
                                             float(line[column_indexes['RecordingTime']])],
                                  [float(prev[column_indexes['GazeRightX']]),
                                   float(line[column_indexes['GazeRightX']])]))
                    # Index 17: Gaze Vector Right Y
                    new_line[column_indexes['GazeRightY']] = str(
                        np.interp(start_ms, [float(prev[column_indexes['RecordingTime']]),
                                             float(line[column_indexes['RecordingTime']])],
                                  [float(prev[column_indexes['GazeRightY']]),
                                   float(line[column_indexes['GazeRightY']])]))
                    # Index 18: Gaze Vector Right Z
                    new_line[column_indexes['GazeRightZ']] = str(
                        np.interp(start_ms, [float(prev[column_indexes['RecordingTime']]),
                                             float(line[column_indexes['RecordingTime']])],
                                  [float(prev[column_indexes['GazeRightZ']]),
                                   float(line[column_indexes['GazeRightZ']])]))

                    # Comment out later -- was just to see which entries in the new file were interpolated
                    # newLine.append('Interpolated')

                    out.write("\t".join(new_line) + "\n")

                    # Update startLine, prev2, and prev references for the next iterations of the for-loop
                    start_line = new_line
                    if last_data:
                        prev2 = prev
                    prev = new_line

                    # This is where the while-loop condition comes into place;
                    # if the difference between the new interpolated entry's RecordingTime and the current line's RecordingTime
                    # is STILL larger than the specified frame-rate interval, then we have to interpolate more data-entries
                    # to fill in this 'large gap' of milliseconds
                    if (float(line[column_indexes['RecordingTime']]) - start_ms) < mspf:
                        small_gap = True
                        if last_data:
                            prev2 = prev
                        prev = line
                        if not delete:
                            out.write("\t".join(line) + "\n")

                            # Just to see how long the whole script takes
                            # print("Finished in {0} seconds".format(time.time() - begin))


def main():
    parser = argparse.ArgumentParser(description='Interpolate data for eye tracker.')
    parser.add_argument('-i', '--input', metavar='/input/file.tsv', type=str,
                        help='path of the input file', required=True)
    parser.add_argument('-o', '--output', metavar='/output/file.tsv', type=str,
                        help='path of the output file', required=True)
    parser.add_argument('--fps', type=int, help='frames per second of interpolation', required=True)
    parser.add_argument('--delete', type=ast.literal_eval, metavar='True|False',
                        help='delete existing values outside of interpolation, defaults to True',
                        required=False, default=True)
    parser.add_argument('--last', type=ast.literal_eval, metavar='True|False',
                        help='keep last data, defaults to True', required=False, default=True)

    args = parser.parse_args()
    interpolate(args.input, args.output, args.fps, delete=args.delete, last_data=args.last)


if __name__ == '__main__':
    main()
