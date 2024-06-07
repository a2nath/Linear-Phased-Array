import json
import argparse
from pathlib import Path, PurePath

# pipe timeout for subprocess
global_timeout    = 60

def findarg(args, key: str) -> bool:
	return key in args and getattr(args, key)

def get_fullpath(output_dir, output_file) -> tuple[Path, str]:
	output_filepath = Path(output_file).resolve()

	if output_filepath.parent != Path(output_dir).resolve():
		output_dir = output_filepath.parent
		output_file = Path(output_file).name

	return Path(output_dir, output_file), output_dir

def convert_lut(input_list: list) -> str:
	return f"\"{' '.join([','.join(map(str, item)) for item in input_list])}\"";

def reparse_binding(input_list: list) -> str:
	return f"\"{' '.join([','.join([':'.join(map(str, timeslot)) for timeslot in timeslots]) for timeslots in input_list])}\"";

# return jason data
def json_data(filename: Path):
	try:
		fp = open(str(filename))
		data = json.load(fp)
		return data, 0
	except ValueError as e:
		return "{}", -1

# validation functoin for the argparser
def file_exists(filepath):
	return os.path.isfile(filepath)

# validation functoin for the argparser
def is_dir(dirpath):
	if os.path.isdir(dirpath):
		return dirpath
	else:
		raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

# validation function for the argparser
def restricted_float(fx):
	try:
		fx = float(fx)
	except ValueError:
		raise argparse.ArgumentTypeError("%r not a floating-point literal" % (fx,))

	if fx < 0.0:
		raise argparse.ArgumentTypeError("%r is not positive"%(fx,))
	return fx
