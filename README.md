# Eye-Track
Eye-tracking device output processing and analysis library.

To run, install the dependencies through ```$ pip install -r requirements.txt``` first. This script needs Python 3.5 to
run. You can view the usage details through ```$ python eyetrack/interpolate.py -h```.

```
usage: interpolate.py [-h] -i /input/file.tsv -o /output/file.tsv --fps FPS
                      [--delete True|False] [--last True|False]

Interpolate data for eye tracker.

optional arguments:
  -h, --help            show this help message and exit
  -i /input/file.tsv, --input /input/file.tsv
                        path of the input file
  -o /output/file.tsv, --output /output/file.tsv
                        path of the output file
  --fps FPS             frames per second of interpolation
  --delete True|False   delete existing values outside of interpolation,
                        defaults to True
  --last True|False     keep last data, defaults to True
```
