The analysis is a single script flybar_process.py written by Charlotte Weaver. It is written in python and requires python 2.7.X and numpy installed.
It relies on the following directory structure: "<top level>/YYYYMMDD/Experiment Name/tracker_files". You can find an example of this structure under trunk/test_data/data. To analyze the test data run

	python flybar_process.py ../test_data/data/
 

It outputs a raw and analyzed file for both the end chamber and the tunnel for each experiment, data from the start chamber is discarded.

	usage: flybar_process.py [-h] [-v] [-d] [--overwrite] directory
	positional arguments:
	  directory      Directory where tracking files are located
	optional arguments:
	  -h, --help     show this help message and exit
	  -v, --verbose  Flag: Increase output verbosity
	  -d, --debug    Flag: Debugging option (also turns on verbosity)
	  --overwrite    Reanalyze all experiments and overwrite previous data
			 (default is to only analyze new data)
