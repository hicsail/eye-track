import numpy as np
import csv
from datetime import datetime, timedelta
# import time
import argparse
import ast
import os
import sys

# Enable relative import for module
PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from eyetrack import config


def interpolate(input_file, output_file, fps, delete=True, file_type='csv'):
    # begin = time.time()
    if fps <= 0:
        print("Cannot have a non-positive value for fps")
        return

    # Turn fps into milliseconds-per-frame
    mspf = (1.0 / fps) * 1000

    with open(input_file, 'r', encoding='UTF-8') as tsv, open(output_file, 'w', encoding='UTF-8') as out:
        if file_type == 'csv':
            read = csv.reader(tsv, delimiter=',')
        else:
            read = csv.reader(tsv, delimiter='\t')

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
        start_ms = 0
        prev = []
        # This was just to parse through the first 39 or so lines of the file for testing purposes, rather than
        #   the entire file
        # i = 0
        for line in read:
            # boolean condition to come into place for the last 'else'-condition
            small_gap = False
            while not small_gap:
                # At the start of a new trial line, demarcated by the 'Separator' in the 'Category'
                # First iteration is skipped, since startMS is initialized to 0 outside the loop
                if line[column_indexes['Category']] == 'Separator':
                    # Write the 'Separator' line, marking the start of the trial, to the new file
                    # -->Sudha said she didn't want the 'Separator' line in the final output file
                    # out.write("\t".join(line)+"\n")

                    # Initialize the prev and startLine references to the first entry after the 'Separator' line,
                    # as well as the startMS for the beginning RecordingTime of the trial
                    prev = next(read)
                    start_ms = float(prev[column_indexes['RecordingTime']])
                    out.write("\t".join(prev) + "\n")
                    small_gap = True
                # If the current read-line's RecordingTime - the startMS of the current time interval < mspf:
                elif (float(line[column_indexes['RecordingTime']]) - start_ms) < mspf:
                    # Update prev reference to the current line, for the next iteration of the for-loop
                    prev = line
                    # If we ARE writing the old data-entries, non-conforming with the constant frame-rate, write them
                    if not delete:
                        out.write("\t".join(line) + "\n")
                    small_gap = True
                # If the current read-line's RecordingTime - the startMS of the current time interval == mspf,
                # AKA, the current data-entry conforms with the specified constant frame-rate -- write the current line
                elif (float(line[column_indexes['RecordingTime']]) - start_ms) == mspf:
                    # The new time interval begins with the RecordingTime of this data-entry
                    start_ms = float(line[column_indexes['RecordingTime']])
                    prev = line
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

                    # For the non-interpolated categories, we have the option of using the previous line's filled-in data
                    # or the current line's filled-in data -- whichever data-entry is closer in time to
                    # the interpolated entry's RecordingTime gets precedence -- and if they're equally distant in time,
                    # use the current line
                    if (start_ms - float(prev[column_indexes['RecordingTime']])) < (
                                float(line[column_indexes['RecordingTime']]) - start_ms):
                        new_line = list(prev)
                    else:
                        new_line = list(line)

                    # Ensure all input values are numeric, otherwise set to 0
                    for i in column_indexes:
                        if i == 'Category' or i == 'TimeOfDay':
                            continue
                        try:
                            float(line[column_indexes[i]])
                        except ValueError:
                            line[column_indexes[i]] = '0'
                        try:
                            float(prev[column_indexes[i]])
                        except ValueError:
                            prev[column_indexes[i]] = '0'

                    # Modify the following index-categories
                    new_line[column_indexes['RecordingTime']] = str(round(start_ms, 1))
                    new_line[column_indexes['TimeOfDay']] = time_of_day
                    # Index 12: Pupil Diameter Right
                    new_line[column_indexes['PupilRight']] = \
                        str(np.interp(start_ms,
                                      [float(prev[column_indexes['RecordingTime']]),
                                       float(line[column_indexes['RecordingTime']])],
                                      [float(prev[column_indexes['PupilRight']]),
                                       float(line[column_indexes['PupilRight']])]))
                    # Index 13: Point of Regard Right X
                    new_line[column_indexes['RegardRightX']] = \
                        str(np.interp(start_ms,
                                      [float(prev[column_indexes['RecordingTime']]),
                                       float(line[column_indexes['RecordingTime']])],
                                      [float(prev[column_indexes['RegardRightX']]),
                                       float(line[column_indexes['RegardRightX']])]))
                    # Index 14: Point of Regard Right Y
                    new_line[column_indexes['RegardRightY']] = \
                        str(np.interp(start_ms,
                                      [float(prev[column_indexes['RecordingTime']]),
                                       float(line[column_indexes['RecordingTime']])],
                                      [float(prev[column_indexes['RegardRightY']]),
                                       float(line[column_indexes['RegardRightY']])]))
                    # Index 16: Gaze Vector Right X
                    new_line[column_indexes['GazeRightX']] = \
                        str(np.interp(start_ms,
                                      [float(prev[column_indexes['RecordingTime']]),
                                       float(line[column_indexes['RecordingTime']])],
                                      [float(prev[column_indexes['GazeRightX']]),
                                       float(line[column_indexes['GazeRightX']])]))
                    # Index 17: Gaze Vector Right Y
                    new_line[column_indexes['GazeRightY']] = \
                        str(np.interp(start_ms,
                                      [float(prev[column_indexes['RecordingTime']]),
                                       float(line[column_indexes['RecordingTime']])],
                                      [float(prev[column_indexes['GazeRightY']]),
                                       float(line[column_indexes['GazeRightY']])]))
                    # Index 18: Gaze Vector Right Z
                    new_line[column_indexes['GazeRightZ']] = \
                        str(np.interp(start_ms,
                                      [float(prev[column_indexes['RecordingTime']]),
                                       float(line[column_indexes['RecordingTime']])],
                                      [float(prev[column_indexes['GazeRightZ']]),
                                       float(line[column_indexes['GazeRightZ']])]))

                    # Comment out later -- was just to see which entries in the new file were interpolated
                    # newLine.append('Interpolated')

                    out.write("\t".join(new_line) + "\n")

                    # Update prev reference for the next iterations of the for-loop
                    prev = new_line

                    # This is where the while-loop condition comes into place;
                    # if the difference between the new interpolated entry's RecordingTime and the current line's RecordingTime
                    # is STILL larger than the specified frame-rate interval, then we have to interpolate more data-entries
                    # to fill in this 'large gap' of milliseconds
                    if (float(line[column_indexes['RecordingTime']]) - start_ms) < mspf:
                        small_gap = True
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
    parser.add_argument('--fps', type=int, help='intended sample rate per seconds', required=True)
    parser.add_argument('--delete', type=ast.literal_eval, metavar='True|False',
                        help='delete existing values outside of interpolation, defaults to True',
                        required=False, default=True)
    parser.add_argument('-t', '--type', type=str, metavar='csv|tsv',
                        help='interpret input file as csv or tsv, defaults to csv',
                        required=False, default='csv')

    args = parser.parse_args()
    interpolate(args.input, args.output, args.fps, delete=args.delete, file_type=args.type)


if __name__ == '__main__':
    main()
