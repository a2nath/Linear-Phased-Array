#!/usr/bin/env python3
import os
import QLearningAgent
from utils.mdp import *


model_bin    = "../Solution/bin/model.exe";
compile_file =  "compile.sh";
limit_sinr   = 24.0
num_episodes = 10;
epsilon      = 0.85;
num_steps    = 10000;
recompile    = False;

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
	parser.add_argument("-b", "--binary", help="Binary to execute for the modelling and simulation", default=model_bin)
	parser.add_argument("-o", "--output_dir",  help="Ouput directory")
	parser.add_argument("-l", "--limit_sinr", \
		help="Lower SINR limit for any handset when calculating reward system for the RL network", default=limit_sinr)
	parser.add_argument("-e", "--num_episodes", help="Number of episodes in the RL simulation", default=num_episodes)
	parser.add_argument("--epsilon", help="Epsilon for greedy explorative agent to check new values", default=epsilon)
	parser.add_argument("-s", "--num_steps", help="Number of steps in each episodes taken by an agent", default=num_steps)
	parser.add_argument("--quiet", help="Debug print off", action='store_true')

	args = parser.parse_args()
	json_datalist = args_init(args);

	# run each scenario file
	for data in json_datalist:

		# set the log Path for logging purposes
		output_file = Path(args.output_dir, "log_" + data[0].name + "." + time.strftime("%Y%m%d-%H%M%S") + ".txt")

		logging.basicConfig(filename=output_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
		console_handler = logging.StreamHandler()
		console_handler.setLevel(logging.INFO)

		logging.getLogger().addHandler(console_handler)
		args.ptx_episodes = [] # across timeslots

		num_timeslots = len(data[1]['base_station_power_dBm_lut']) # 1

		env = MDPEnv(args.binary, args.limit_sinr, not args.quiet, logging.getLogger(). data[1])

		# get powers from previous timeslot

		agent = QLearningAgent(env.action_space, env.observation_space.shape[0])

		for episode in range(args.num_episodes):
			start_time = timeit.default_timer()

			# erase the state (as opposed to advance())
			state = env.reset(args.timeslot, data)
			total_reward = [0] * num_timeslots

			for timeslot in num_timeslots:

				done = False
				steps = 0
				total_reward[timeslot] = 0

				while not done and steps < args.num_steps:
					action = agent.choose_action(state)
					next_state, reward, done, _ = env.step(action)
					agent.learn(state, action, reward, next_state)
					state = next_state
					total_reward[timeslot] += reward
					steps += 1

				logging.info(f"state {env.results}, reward {total_reward[timeslot]}")

				# critical for simulating more than 1 timeslot
				env.advance_simulation()

			agent.update_epsilon()
			print(f"Episode {episode + 1}: Total Reward: {total_reward}")

			end_time = timeit.default_timer()


		print("-------------------------------------------------------")
		print(f"{data[0].name} test done in {end_time - start_time}")

if __name__ == '__main__':
	main()