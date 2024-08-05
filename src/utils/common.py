import json
import argparse
from typing import Union, TextIO
from pathlib import Path, PurePath
import logging

# pipe timeout for subprocess
global_timeout = 60
global_c_speed = 3e8

# give file handle instead of file to avoid opening and closing
def log_print(log_handle: Union[TextIO, str], text: str, mode: str = 'a+'):
	if isinstance(log_handle, str):
		with open(log_handle, mode) as file:
			file.write(text + '\n')
	else:
		log_handle.write(text + '\n')
	log_handle.flush()

def log_printline(log: TextIOWrapper, text: str, mode = 'a+'):
	log_print(log, text + '\n', mode)

def log_printfor(log: TextIOWrapper, array, mode = 'a+', delim = ','):
	newstr = array[0]
	for item in array[1:]:
		newstr += delim + item
	log_printline(log, newstr, mode)

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

def reparse_binding(input_list: [[[]]]) -> str:
	return ' '.join([','.join([':'.join([str(node) for node in timeslot]) for timeslot in timeslots]) for timeslots in input_list])


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



class Argument:
	# val        : value of the argument
	# mutable    : can be changed between steps RL
	# parameters : counts as 1 or more parameter(s)

	@property
	def value(self):
		return self._val

	@value.setter
	def value(self, val):
		self._val = val


	@property
	def len(self):
		return self._num_parameters

	@property
	def mutable(self):
		return self._mutable

	def __init__(self, val, mutable = False, parameters = 1):
		self._val = val
		self._mutable = mutable
		self._num_parameters = parameters
