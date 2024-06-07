#!/usr/bin/env python3
import os
import QLearningAgent
from utils.mdp import *

model_bin          = "../Solution/bin/model.exe";
compile_file       =  "compile.sh";
sinr_limit         = 24.0
number_of_episodes = 1;
epsilon            = 0.1;
number_of_steps    = 1000;
recompile          = False;

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
	parser.add_argument("-l", "--snr_limit", \
		help="Lower SINR limit for any handset when calculating reward system for the RL network", default=sinr_limit)
	parser.add_argument("-e", "--episode", help="Number of episodes in the RL simulation", default=number_of_episodes)
	parser.add_argument("--epsilon", help="Epsilon for greedy explorative agent to check new values", default=epsilon)
	parser.add_argument("-s", "--step", help="Number of steps in each episodes taken by an agent", default=number_of_steps)
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