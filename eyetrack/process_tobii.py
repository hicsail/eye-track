import argparse


# Slight implementation note: I actually accidentally wrote all the string-values to the outputfile.csv as 'some_string', whereas
# I didn't notice after the fact that all such string values in the original input files are written as '"some_string"', and written to
# the outputfile.csv that way...
def process_tobii(tobiifile, aoifile, phasefile, outputfile, condition, trialorder):
    # Read the files in as 2D arrays
    with open(tobiifile) as t, open(aoifile) as a, open(phasefile) as phaseTiming:
        tobii = [line.split("\t")[:-1] for line in t]
        aoi = [line.strip().split("\t") for line in a]
        phase = [line.strip().split("\t") for line in phaseTiming]

    # I declared a ton of variables that would make it easier (in case certain column-indices switched around in different file formats)
    # to simply store the given index of a specific column at the beginning (from the column-headers), and then continually use that stored index throughout,
    # by way of just a key into a dict
    tobii_indices = {}
    aoi_indices = {}
    phase_indices = {}

    # Create the additional headers for the created-columns now for the tobii-file
    tobii[0] += ["X", "Condition", "TrialOrder", "Trial", "Phase", "FramesFromMovieOnset", "FramesFromPhaseOnset",
                 "TimeFromMovieOnset", "TimeFromPhaseOnset",
                 "Subphase", "FramesFromSubphaseOnset", "TimeFromSubphaseOnset", "GazeXAvg", "GazeYAvg", "SceneName",
                 "Position", "SceneType"]
    for i in range(len(tobii[0])):
        # rename some columns
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
        # store the index of that column
        tobii_indices[tobii[0][i]] = i
    # refer to them later easily by memorable variable-names
    movie = tobii_indices["MovieName"]
    tobii_participant = tobii_indices["ParticipantName"]
    event = tobii_indices["StudioEvent"]
    gaze_right_x = tobii_indices["GazePointRightX"]
    gaze_right_y = tobii_indices["GazePointRightY"]
    gaze_left_x = tobii_indices["GazePointLeftX"]
    gaze_left_y = tobii_indices["GazePointLeftY"]
    validity_left = tobii_indices["ValidityLeft"]
    validity_right = tobii_indices["ValidityRight"]
    # cond = tobii_indices["Condition"]
    # t_order = tobii_indices["TrialOrder"]
    tobii_trial = tobii_indices["Trial"]
    tobii_phase = tobii_indices["Phase"]
    f_movie = tobii_indices["FramesFromMovieOnset"]
    f_phase = tobii_indices["FramesFromPhaseOnset"]
    t_movie = tobii_indices["TimeFromMovieOnset"]
    t_phase = tobii_indices["TimeFromPhaseOnset"]
    tobii_subphase = tobii_indices["Subphase"]
    f_subphase = tobii_indices["FramesFromSubphaseOnset"]
    t_subphase = tobii_indices["TimeFromSubphaseOnset"]
    gaze_x_avg = tobii_indices["GazeXAvg"]
    gaze_y_avg = tobii_indices["GazeYAvg"]
    tobii_scene = tobii_indices["SceneName"]
    tobii_pos = tobii_indices["Position"]
    tobii_type = tobii_indices["SceneType"]

    # same for aoi-file
    for i in range(len(aoi[0])):
        aoi_indices[aoi[0][i]] = i
    aoi_trial = aoi_indices["Trial"]
    scene = aoi_indices["SceneName"]
    pos = aoi_indices["Position"]
    aoi_type = aoi_indices["SceneType"]
    max_aoi = aoi_indices["Top Right (x,y)"]
    min_aoi = aoi_indices["Bottom Left (x,y)"]

    # This was to store the entries of the aoi-file by using the filename as the key, and the row's entry
    # as the value -- such that I wouldn't have to loop over the entire aoi-file for every row in the tobii-file
    aoi_trials = {}
    for i in range(1, len(aoi)):  # skip the first line of column-headers
        trialkey = aoi[i][aoi_trial]
        if trialkey in aoi_trials:
            aoi_trials[trialkey] += [aoi[i]]
        else:
            aoi_trials[trialkey] = [aoi[i]]

    for i in range(len(phase[0])):
        phase_indices[phase[0][i]] = i
    filename = phase_indices['Filename']
    trial = phase_indices['Trial']
    phases = phase_indices['Phase']
    subphase = phase_indices['Subphase']
    start = phase_indices['Starting_time']

    # Same concept here with the phasefile as with the aoi-file, except I also store the (Subphase, Starting_time) as the value,
    # with the key being the 'Trial' category
    phase_filenames = {}
    subphase_times = {}
    for i in range(1, len(phase)):  # skip the first line of column-headers
        filekey = phase[i][filename]
        sublist = phase[i]
        to_add = [sublist[trial], sublist[phases]]
        if filekey in phase_filenames:
            phase_filenames[filekey] += [to_add]
        else:
            phase_filenames[filekey] = [to_add]

        trialkey = phase[i][trial]
        if trialkey in subphase_times:
            subphase_times[trialkey] += [(phase[i][subphase], phase[i][start])]
        else:
            subphase_times[trialkey] = [(phase[i][subphase], phase[i][start])]

    hierarchy = {}
    index = 1  # skip the first line of column-headers
    while index != len(tobii):
        # Essentially a while-loop that iterates over the tobii-file, but the incremented indexes are kept consistent
        # as I continually pop any rows that have the following condition:
        if tobii[index][event] == "MovieStart" or tobii[index][event] == "MovieEnd":
            tobii.pop(index)
        else:
            # to_append is the enhanced row-entry with the extra columns, all with default '' values for now
            to_append = tobii[index] + ['', condition, trialorder] + [''] * 14

            # Simultaneously fill in the Trial and Phase categories for each row, while filling in any leading NA-gaps
            cur_movie = tobii[index][movie]
            if cur_movie in phase_filenames:
                # 0 for middle-index because I think, for a given filename, the trial and phase are consistent throughout the rows
                to_append[tobii_trial] = phase_filenames[cur_movie][0][0]
                to_append[tobii_phase] = phase_filenames[cur_movie][0][1]
            elif index == 1:
                pass
            else:
                to_append[tobii_trial] = tobii[index - 1][tobii_trial]

            # Here, what I'm doing is filtering out "blocks" of rows by ParticipantName --> Condition --> MovieName/Trial --> Phase --> Subphase,
            # by simply creating a hierarchy of nested dictionaries that each row-entry assigns itself to, incrementing the counter or becoming
            # the first (so the 0th frame) of that "block" of rows --> This is simpler than going through ParticipantName --> ... --> MovieName/Trial,
            # ParticipantName --> ... --> Subphase, and calculating separately the 'FramesFromSubphaseOnset' in several different for-loops
            cur_partic = to_append[tobii_participant]
            cur_trial = to_append[tobii_trial]
            cur_phase = to_append[tobii_phase]
            if cur_partic != '' and condition != '' and cur_movie != '':
                if cur_partic not in hierarchy:
                    hierarchy[cur_partic] = {condition: {cur_movie: 0}}
                else:
                    if condition not in hierarchy[cur_partic]:
                        hierarchy[cur_partic][condition] = {cur_movie: 0}
                    else:
                        if cur_movie not in hierarchy[cur_partic][condition]:
                            hierarchy[cur_partic][condition][cur_movie] = 0
                        else:
                            hierarchy[cur_partic][condition][cur_movie] += 1
                # Write to the 'FramesFromMovieOnset' and 'Time...' categories
                to_append[f_movie] = hierarchy[cur_partic][condition][cur_movie]
                to_append[t_movie] = round(float(to_append[f_movie]) * 1000.0 / 60.0)  # it just says "The data is collected at 60Hz"

                if cur_trial != '' and cur_phase != '':
                    if cur_trial not in hierarchy[cur_partic][condition]:
                        hierarchy[cur_partic][condition][cur_trial] = {cur_phase: {'framesFromPhaseOnset': 0}}
                    else:
                        if cur_phase not in hierarchy[cur_partic][condition][cur_trial]:
                            hierarchy[cur_partic][condition][cur_trial][cur_phase] = {'framesFromPhaseOnset': 0}
                        else:
                            hierarchy[cur_partic][condition][cur_trial][cur_phase]['framesFromPhaseOnset'] += 1
                    # Write to the 'FramesFromPhaseOnset' and 'Time...' categories
                    to_append[f_phase] = hierarchy[cur_partic][condition][cur_trial][cur_phase]['framesFromPhaseOnset']
                    to_append[t_phase] = round(
                        float(to_append[f_phase]) * 1000.0 / 60.0)  # it just says "The data is collected at 60Hz"

                    if cur_trial in subphase_times:
                        # I iterate, for each row with a Subphase, through the a subset of entries in the phasefile to
                        # see the closest Starting_time it comes after
                        cur_timestamp = float(to_append[t_phase])
                        matched_subphases = subphase_times[cur_trial]
                        closest_time = 0
                        for x in range(1, len(matched_subphases)):
                            closest_diff = cur_timestamp - float(matched_subphases[closest_time][1]) * 1000.0
                            cur_diff = cur_timestamp - float(matched_subphases[x][1]) * 1000.0
                            if 0.0 <= cur_diff <= closest_diff:
                                closest_time = x
                        to_append[tobii_subphase] = matched_subphases[closest_time][0]

                        cur_subphase = to_append[tobii_subphase]
                        if cur_subphase not in hierarchy[cur_partic][condition][cur_trial][cur_phase]:
                            hierarchy[cur_partic][condition][cur_trial][cur_phase][cur_subphase] = 0
                        else:
                            hierarchy[cur_partic][condition][cur_trial][cur_phase][cur_subphase] += 1
                        # Write to the 'FramesFromSubphaseOnset' and 'Time...' categories
                        to_append[f_subphase] = hierarchy[cur_partic][condition][cur_trial][cur_phase][cur_subphase]
                        to_append[t_subphase] = round(
                            float(to_append[f_subphase]) * 1000.0 / 60.0)  # it just says "The data is collected at 60Hz"
                else:
                    # These are just to stay consistent with how the original R-script wrote 0 as the default instead of NA, no matter what
                    to_append[f_subphase] = 0
                    to_append[t_subphase] = 0
            else:
                to_append[f_movie] = 0
                to_append[t_movie] = 0
                to_append[f_phase] = 0
                to_append[t_phase] = 0
                to_append[f_subphase] = 0
                to_append[t_subphase] = 0

            # Here, we compute the GazeX/YAvg
            c_validity_left = float(to_append[validity_left])
            c_validity_right = float(to_append[validity_right])
            # Some rows only had the right-eye GazePoint coordinates, some only had the left, some had neither -- so we only convert to float if they exist
            if to_append[gaze_right_x] != '' and to_append[gaze_right_y] != '':
                c_gaze_right_x = float(to_append[gaze_right_x])
                c_gaze_right_y = float(to_append[gaze_right_y])
            if to_append[gaze_left_x] != '' and to_append[gaze_left_y] != '':
                c_gaze_left_x = float(to_append[gaze_left_x])
                c_gaze_left_y = float(to_append[gaze_left_y])
            # Check validity and assign/compute GazeX/YAvg
            if (c_validity_left > 1 or c_validity_left < 0) and (c_validity_right > 1 or c_validity_right < 0):
                to_append[gaze_x_avg] = -999
                to_append[gaze_y_avg] = -999
            elif (c_validity_left > 1 or c_validity_left < 0) and (1 >= c_validity_right >= 0):
                # take the right side
                to_append[gaze_x_avg] = c_gaze_right_x
                to_append[gaze_y_avg] = c_gaze_right_y
            elif (c_validity_right > 1 or c_validity_right < 0) and (1 >= c_validity_left >= 0):
                # take the left side
                to_append[gaze_x_avg] = c_gaze_left_x
                to_append[gaze_y_avg] = c_gaze_left_y
            else:
                # avg both sides
                to_append[gaze_x_avg] = (c_gaze_left_x + c_gaze_right_x) / 2.0
                to_append[gaze_y_avg] = (c_gaze_left_y + c_gaze_right_y) / 2.0
            # We only attempt to write the SceneName, Position, and SceneType the following under these conditions:
            if to_append[gaze_x_avg] != -999 and cur_trial != '' and cur_trial in aoi_trials:
                matched_aoi = aoi_trials[cur_trial]
                for area in matched_aoi:
                    # Quite annoying -- parsing (x,y) coordinates where the numbers themselves can be written with commas if they're > 1000 -- like 45 lines of code for just this
                    max_coords = area[max_aoi].strip()
                    min_coords = area[min_aoi].strip()
                    max_commas = max_coords.count(',')
                    if max_commas == 1:
                        max_x_y = max_coords.split(',')
                        x_max = int(max_x_y[0])
                        y_max = int(max_x_y[1])
                    elif max_commas == 2:
                        max_x_y = max_coords.split(',')
                        temp1 = int(max_x_y[0]) * (10 ** len(max_x_y[1])) + int(max_x_y[1])
                        temp2 = int(max_x_y[1]) * (10 ** len(max_x_y[2])) + int(max_x_y[2])
                        if (temp1 - int(max_x_y[2])) > (temp2 - int(max_x_y[0])):
                            # use temp2
                            x_max = int(max_x_y[0])
                            y_max = temp2
                        else:
                            # use temp1
                            x_max = temp1
                            y_max = int(max_x_y[2])
                    else:  # (For the case of 3 commas)
                        max_x_y = max_coords.split(',')
                        x_max = int(max_x_y[0]) * (10 ** len(max_x_y[1])) + int(max_x_y[1])
                        y_max = int(max_x_y[2]) * (10 ** len(max_x_y[3])) + int(max_x_y[3])

                    min_commas = min_coords.count(',')
                    if min_commas == 1:
                        min_x_y = min_coords.split(',')
                        x_min = int(min_x_y[0])
                        y_min = int(min_x_y[1])
                    elif min_commas == 2:
                        min_x_y = min_coords.split(',')
                        temp1 = int(min_x_y[0]) * (10 ** len(min_x_y[1])) + int(min_x_y[1])
                        temp2 = int(min_x_y[1]) * (10 ** len(min_x_y[2])) + int(min_x_y[2])
                        if (temp1 - int(min_x_y[2])) > (temp2 - int(min_x_y[0])):
                            # use temp2
                            x_min = int(min_x_y[0])
                            y_min = temp2
                        else:
                            # use temp1
                            x_min = temp1
                            y_min = int(min_x_y[2])
                    else:  # (For the case of 3 commas)
                        min_x_y = min_coords.split(',')
                        x_min = int(min_x_y[0]) * (10 ** len(min_x_y[1])) + int(min_x_y[1])
                        y_min = int(min_x_y[2]) * (10 ** len(min_x_y[3])) + int(min_x_y[3])
                    # Since I don't know the format around the pixel-coordinates of the corners of the area-of-interest on the eye-tracker screen, apparently
                    # the "Bottom-Left" y-pixel-coordinate value can be greater than the "Top-Right" y-pixel-coordinate value
                    if max(x_max, x_min) >= to_append[gaze_x_avg] >= min(x_max, x_min) \
                            and max(y_max, y_min) >= to_append[gaze_y_avg] >= min(y_max, y_min):
                        to_append[tobii_scene] = area[scene]
                        to_append[tobii_pos] = area[pos]
                        to_append[tobii_type] = area[aoi_type]
            # Otherwise, only in the case where ValidityLeft and ValidityRight were unacceptable, do we assign the value 'TrackLost''
            elif to_append[gaze_x_avg] == -999:
                to_append[tobii_scene] = "TrackLoss"
                to_append[tobii_pos] = "TrackLoss"
                to_append[tobii_type] = "TrackLoss"

            tobii[index] = to_append
            index += 1
    # For some reason, only the 'StudioEvent' column could be left with '' values; all other '' turned into NA
    for i in range(len(tobii)):
        for j in range(len(tobii[i])):
            if j == event:
                continue
            elif tobii[i][j] == '':
                tobii[i][j] = 'NA'
            elif type(tobii[i][j]) != str:
                tobii[i][j] = str(tobii[i][j])  # For numbers --> strings
    with open(outputfile, 'w') as out:
        for line in tobii:
            out.write(",".join(line) + "\n")  # Comma-separated-values


# Checks whether a string represents a number
def is_numeric(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def check_python_process_tobii(r_outputfile, python_outputfile):
    def print_differences():
        print("Not equal at cell (1-indexed): ({0}, {1})".format(i + 1, column_char))
        print("R-Output: {0}".format(r_file[i][j]))
        print("Python-Output: {0}".format(py_file[i][j]))

    with open(r_outputfile) as Rout, open(python_outputfile) as Pout:
        r_file = [line.strip().split(",") for line in Rout]
        py_file = [line.strip().split(",") for line in Pout]
    if len(py_file) != len(r_file):
        print("Lengths aren't equal")
    else:
        for i in range(len(r_file)):
            for j in range(len(r_file[i])):
                if j < 26:
                    # column_char is supposed to represent the first 52 column-headers Microsoft Excel displays for a .csv file (so A, B, C, ... AA, AB, AC, ...)
                    column_char = str(chr(j + 97)).upper()
                else:
                    column_char = "A" + str(chr((j % 26) + 97)).upper()
                if i == 0 and j == 3:  # This was the cell for 'Gender' -- somewhat strange how the original R-script changed the column-header
                    continue
                # If the value in the column is a string and not NA
                elif not is_numeric(py_file[i][j]) and py_file[i][j] != "NA":
                    # Put back in the '"some_string"' double-quotes
                    if r_file[i][j] != py_file[i][j]:
                        print_differences()
                        return
                elif is_numeric(py_file[i][j]):
                    # Cast both values in the R-output and the Python-output to float values, for precision purposes (as '4.10' != '4.1')
                    if float(r_file[i][j]) != float(py_file[i][j]):
                        print_differences()
                        return
                elif r_file[i][j] != py_file[i][j]:
                    print_differences()
                    return
    print("Files are identical")


def main():
    parser = argparse.ArgumentParser(description='Parse data for eye tracker.')
    parser.add_argument('-i', '--input', metavar='/input/file.tsv', type=str,
                        help='path of the input file', required=True)
    parser.add_argument('-a', '--aoi', metavar='/input/aiofile.tsv', type=str,
                        help='path of the input AOI file', required=True)
    parser.add_argument('-p', '--phase', metavar='/input/phasefile.tsv', type=str,
                        help='path of the input Phase file', required=True)
    parser.add_argument('-o', '--output', metavar='/output/file.tsv', type=str,
                        help='path of the output file', required=True)
    parser.add_argument('--condition', metavar='Intransitive', type=str,
                        help='condition', required=False, default='Intransitive')
    parser.add_argument('--order', metavar='Backward', type=str,
                        help='order', required=False, default='Backward')

    args = parser.parse_args()
    process_tobii(args.input, args.aoi, args.phase, args.output, args.condition, args.order)


if __name__ == '__main__':
    main()
