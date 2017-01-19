# Eye-Track
Eye-tracking device output processing and analysis library.

To run, install the dependencies through ```$ pip install -r requirements.txt``` first. This script needs Python 3.5 to
run.
Two tools are available: ```process_tobii.py``` and ```interpolate.py```. Use process_tobii.py to merge eye tracker
data with the AIO and Phase files. Use interpolate.py to obtain a constant sample rate.

You can view the usage details through ```$ python eyetrack/interpolate.py -h```.
```
usage: interpolate.py [-h] -i /input/file.tsv -o /output/file.tsv --fps FPS
                      [--delete True|False] [-t csv|tsv]

Interpolate data for eye tracker.

optional arguments:
  -h, --help            show this help message and exit
  -i /input/file.tsv, --input /input/file.tsv
                        path of the input file
  -o /output/file.tsv, --output /output/file.tsv
                        path of the output file
  --fps FPS             intended sample rate per seconds
  --delete True|False   delete existing values outside of interpolation,
                        defaults to True
  -t csv|tsv, --type csv|tsv
                        interpret input file as csv or tsv, defaults to csv
```

You can view the usage details through ```$ python eyetrack/process_tobii.py -h```.
```
usage: process_tobii.py [-h] -i /input/file.tsv -a /input/aiofile.tsv -p
                        /input/phasefile.tsv -o /output/file.tsv
                        [--condition Intransitive] [--order Backward]
                        [-t csv|tsv]

Parse data for eye tracker.

optional arguments:
  -h, --help            show this help message and exit
  -i /input/file.tsv, --input /input/file.tsv
                        path of the input file
  -a /input/aiofile.tsv, --aoi /input/aiofile.tsv
                        path of the input AOI file
  -p /input/phasefile.tsv, --phase /input/phasefile.tsv
                        path of the input Phase file
  -o /output/file.tsv, --output /output/file.tsv
                        path of the output file
  --condition Intransitive
                        condition
  --order Backward      order
  -t csv|tsv, --type csv|tsv
                        interpret input files as csv or tsv, defaults to csv
```
