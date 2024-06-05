#!/usr/bin/env python3

import os
import sys
import argparse
from subprocess import run
from pathlib import Path, PurePath
import argparse
import json
import subprocess
import time, timeit
import QLearningAgent
import utils.common
import pdb

model_bin = "../Solution/bin/model.exe";
compile_file = "compile.sh";
default_SINR = 24.0
recompile = False;

number_of_episodes = 1000;
global_timeout = 60

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




class MDP:

	def get_output(self, args):
		global global_timeout
		sinr_data = []
		start_time = timeit.default_timer()
		if Path(self.model_bin).is_file():
			print("model path ok")

		print([self.model_bin] + args)

		process = subprocess.Popen([self.model_bin] + args, \
			stdout=subprocess.PIPE, \
			stderr=subprocess.PIPE, bufsize=1)

		pdb.set_trace()
		while True:
			stdout_line_str = process.stdout.readline()
			if not stdout_line_str and process.poll() is not None:
				break

			elif not stdout_line_str:
				while True:
					stdout_line_str = process.stderr.readline()

					if not stdout_line_str and process.poll() is not None:
						break

					# Read and print the output line by line
					if stdout_line_str and self.debug_flag:
						print(stdout_line_str.rstrip(), flush=True)
					break

			try:
				stdout_line_str = stdout_line_str.rstrip()
				sinr_data.append(stdout_line_str.split(' '))
			except UnicodeDecodeError as e:
				# Handle decoding error gracefully
				print(f"Error decoding line from stdout: {e}")
				print(f"Problematic byte sequence: {stdout_line_str}")
				continue
			except KeyboardInterrupt:
				process.kill()
				outs, errs = process.communicate(timeout = global_timeout)
				print("----------------------KILLED---------------------------")
				print(f"Process interrupted.\nouts = {outs.decode('utf-8').rstrip()}\nerrs = {errs.decode('utf-8').rstrip()}")

			if global_timeout is not None and float(timeit.default_timer() - start_time) > float(global_timeout):
				process.kill()
				outs, errs = process.communicate(timeout=global_timeout)
				print("----------------------KILLED---------------------------")
				print(f"Process timed out after {global_timeout} seconds\nouts = {outs.decode('utf-8').rstrip()}\nerrs = {errs.decode('utf-8').rstrip()}")


		return sinr_data

	# input:  jsondata: list, output_logfile: Path
	# output: return_code, [sinr_values]
	def action(self, input_args, output_logfile: Path):
		output = get_output(model_bin, input_args)

		print(sirn_data)

		return output

		#print(
		#	f"Input:\t{input_args}",
		#	f"Output:\t{sirn_data}\n",
		#	file=output_logfile,
		#	flush=True,
		#)

	def init_args(self, jdata: list):

		# validate
		j_antenna_txpower = jdata['base_station_power_dBm_lut']
		j_scan_angle      = jdata['base_station_scan_alpha_deg_lut']
		j_ms_selection    = jdata['base_to_mobile_station_id_selection_lut']

		timeslots = 0

		if len(j_antenna_txpower) == len(j_scan_angle) and len(j_scan_angle) == len(j_ms_selection):
			timeslots = len(j_antenna_txpower)
		else:
			raise ValueError("The number of timeslots is mismatch between power, scan, ms:bs id mapping.")

		self.sim_args = [\
			'--frequency',				str(jdata['frequency']),\
			'--bandwidth',				str(jdata['bandwidth']),\
			'--symbolrate',				str(jdata['symbolrate']),\
			'--blockpersymbol',			str(jdata['blockpersymbol']),\
			'--height',					str(jdata['height']),\
			'--ms_grx',					str(jdata['ms_grx_dB']),\
			'--system_noise',			str(jdata['system_noise']),\
			'--mobile_stations',		str(jdata['mobile_stations']),\
			'--base_stations',			str(jdata['base_stations']),\
			'--timeslots',			    str(timeslots),\
			'--slimit',			        str(jdata['sinr_limit']),\
			'--base_station_theta_c',	','.join(map(str, jdata['base_station_theta_c'])),\
			'--base_station_location',		convert_lut(jdata['base_station_location']),\
			'--mobile_station_location',	convert_lut(jdata['mobile_station_location']),\
			'--base_station_antenna_counts',			','.join(map(str, jdata['base_station_antenna_counts'])),\
			'--base_station_power_range_dBm',			','.join(map(str, jdata['base_station_power_range_dBm'])),\
			'--base_station_scan_angle_range_deg',		','.join(map(str, jdata['base_station_scan_angle_range_deg'])),\
			'--antenna_spacing',		','.join(map(str, jdata['antenna_spacing'])),\
			'--antenna_dims',			','.join(map(str, jdata['antenna_dims'])),\
			# suggested starting point
			'--antenna_txpower',            convert_lut(j_antenna_txpower),\
			'--scan_angle',	                convert_lut(j_scan_angle),\
			'--ms_selection',	            reparse_binding(jdata['base_to_mobile_station_id_selection_lut']),\
		]


		return args

	# state_space_size =  # Define the size based on your state representation
	# action_space_size =  # Define the size based on your action space
	# agent = QLearningAgent(state_space_size, action_space_size)

	# a [step, S] will update the args as per the [rewards, R] system and take [action, A]
	def task(self):
		upgdate_args()
		output = action()

		average_sinr = 0;
		for item in output:
			average_sinr += item.sinr;


		#expected_output =
		#delta  = expected_output - output
		#reward =

		self.Q = average_sinr / len(output)

		#print(sirn_data)




	def launch_scenario(self, data):

		for episode in range(self.episode):
			self.init_args(data)

			for step in range(self.steps):
				task()


	def __init__(self, args, debug = False):

		# validate the executatable
		self.debug_flag       = debug
		self.model_bin        = default_model_bin
		self.episode          = args.episode
		self.steps            = args.step
		self.output_dir       = args.output_dir
		self.reward_limit     = args.snr_limit
		self.timeout          = None
		self.sim_args         = None
		#else:
		#	raise FileNotFoundError(f"Executable file has a problem or does not exist {default_model_bin}")


def args_init(args):
	json_datalist = [] # stores json data, rather than files

	if findarg(args, 'filename'):
		args.filename, output_dir = get_fullpath(args.input_dir, args.filename)

		if not findarg(args, 'output_dir'):
			args.output_dir = output_dir # set the output dir if not set

		if args.filename.is_file(): # skip checking .json in the filename
			data, retval = json_data(args.filename)
			if retval == 0:
				json_datalist.append([args.filename, data])

	elif findarg(args, 'input_dir'):

		if not findarg(args, 'output_dir'):
			args.output_dir = args.input_dir # set the output dir if not set

		for filename in os.listdir(args.input_dir):
			filename = Path(filename)
			if filename.is_file(): # skip checking .json in the filename
				data, retval = json_data(filename)
				if retval == 0:
					json_datalist.append([filename, data])

	if len(json_datalist) == 0:
		print("There were no files to process")
		exit(0)

	print("Found ", len(json_datalist), " files")

	# make the directory if missing
	if args.output_dir != args.input_dir:
		Path(args.output_dir).mkdir(parents=True, exist_ok=True)

	# print out the settings before running
	if not args.quiet:
		print("-----------------------SETTINGS------------------------")

		arguments = vars(args);
		for key in arguments:
			value = getattr(args, key)
			if (key and type(value) != bool and value != None) or type(value) == bool:
				print(key, '\t', value)

		print("-----------------------STARTING------------------------")

	return json_datalist


def main():

	parser = argparse.ArgumentParser("Performs Reinforcement Learning to check for the most optimized location of nodes in the network", add_help=True)
	parser.add_argument("-f", "--filename", help="Name of the input json file")
	parser.add_argument("-i", "--input_dir", help="Input directory where input json files are", default=os.getcwd())
	parser.add_argument("-o", "--output_dir",  help="Ouput directory")
	parser.add_argument("-l", "--snr_limit", help="Lower SINR limit for any handset when calculating reward system for the RL network", default=default_SINR)
	parser.add_argument("-e", "--episode", help="Number of episodes in the RL simulation", default=1)
	parser.add_argument("-s", "--step", help="Number of steps in each episodes taken by an agent", default=1000)
	parser.add_argument("--quiet", help="Debug print off", action='store_true')

	args = parser.parse_args()
	json_datalist = args_init(args);
	mdp_process = MDP(args, not args.quiet)

	# run each scenario file
	for data in json_datalist:
		start_time = timeit.default_timer()

		output_file = str(Path(args.output_dir, data[0].name + "." + time.strftime("%Y%m%d-%H%M%S") + ".txt"))
		mdp_process.launch_scenario(data[1], output_file)

		end_time = timeit.default_timer()
		print("-------------------------------------------------------")
		print(f"{data[0].name} test done in {end_time - start_time}")

if __name__ == '__main__':
	main()