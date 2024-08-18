#!/usr/bin/env python3
import os, sys
import QLearningAgent
from utils.mdp import *

model_bin    = "../Solution/bin/model.exe"
comile_script = "../compile.sh"
compile_script =  "compile.sh"
min_snr = 24.0
num_episodes = 10
learning_rate = 0.1
discount_factor = 0.99
epsilon_exp = 1; # exploration_rate
epsilon_min = 0.01
epsilon_max = 1
epsilon_dec = 0.001
num_steps   = 10000
recompile   = False

def args_init(args):
	# check for binary and recompile if needed
	if file_exists(compile_script) and args.compile:
		exec(f"./{compile_script}")
	elif args.compile:
		sys.exit(f"Compile script not found or not valid{args.compile}")

	# stores json data structure, rather than files
	if file_exists(args.binary):
		json_datalist = []

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

	else:
		sys.exit("Binary doesn't exist. Pass in the --compile flag and rerun the script.")


def main():
	parser = argparse.ArgumentParser("Performs Reinforcement Learning to check for the most optimized location of nodes in the network", add_help=True)
	parser.add_argument("-f", "--filename", help="Name of the input json file")
	parser.add_argument("-i", "--input_dir", help="Input directory where input json files are", default=os.getcwd())
	parser.add_argument("-b", "--binary", help="Binary to execute for the modelling and simulation", default=model_bin)
	parser.add_argument("-o", "--output_dir",  help="Ouput directory")
	parser.add_argument("-s", "--min_snr", \
		help="Lower SINR limit for any handset when calculating reward system for the RL network", default=min_snr)
	parser.add_argument("-p", "--num_episodes", help="Number of episodes in the RL simulation", default=num_episodes)
	parser.add_argument("-d", "--discount_factor", help="The weightage to give to the new adjustment of the q table", default=discount_factor)
	parser.add_argument("-l", "--learning_rate", help="Learning rate or the weightage to give to old values", default=learning_rate)
	parser.add_argument("-e", "--epsilon", help="Epsilon for greedy explorative agent to check new values", default=epsilon_exp)
	parser.add_argument("-c", "--epsilon_decay", help="Epsilon value decay rate for each episode", default=epsilon_dec)
	parser.add_argument("-t", "--num_steps", help="Number of steps in each episodes taken by an agent", default=num_steps)
	parser.add_argument("-n", "--dtype", help=f"ML precision for Q table", default=np.zeros(1).dtype)#{','.join(supported_precisions)}",\
	#	default=np.zeros(1).dtype, choices=supported_precisions)
	parser.add_argument("--compile", help="Recompile the binary and run the script", action='store_true')
	parser.add_argument("--quiet", help="Debug print off", action='store_true')

	args = parser.parse_args()
	json_datalist = args_init(args);
	precision = np.dtype(args.dtype)



	# run each scenario file
	for data in json_datalist:

		# set the log Path for logging purposes
		output_file = Path(args.output_dir, "log_" + data[0].name + "." + time.strftime("%Y%m%d-%H%M%S") + ".txt")

		logging.basicConfig(filename=output_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
		console_handler = logging.StreamHandler()
		console_handler.setLevel(logging.INFO)

		logging.getLogger().addHandler(console_handler)
		states_per_episode = {} # across timeslots
		rewards_per_episodes = []

		num_timeslots = len(data[1]['base_station_power_dBm_lut']) # 1

		env = MDPEnv(args.binary, args.min_snr, not args.quiet, logging.getLogger(), data[1])

		# get powers from previous timeslot
		action_space = spaces.MultiDiscrete([3] * env.num_parameters, dtype=np.uint8)
		agent = QLearningAgent(action_space, env.observation_space, \
						 args.learning_rate, args.discount_factor, args.epsilon, args.epsilon_decay, epsilon_min, epsilon_max, precision)

		for episode in range(args.num_episodes):
			start_time = timeit.default_timer()

			# erase the state (as opposed to advance())
			total_reward = [0] * num_timeslots
			episode_state = []

			for timeslot in range(num_timeslots):

				state = env.reset(timeslot)
				done = False
				steps = 0
				#total_reward[timeslot] = 0

				while not done and steps < args.num_steps:
					action = agent.choose_action(state)
					next_state, reward, done = env.step(action)
					agent.learn(state, action, reward, next_state)
					state = next_state
					total_reward[timeslot] += reward
					steps += 1

				logging.info(f"state {env.results}, reward {total_reward[timeslot]}")

				# critical for simulating more than 1 timeslot
				episode_state.append(env.state())
				env.advance_simulation()

			agent.update_epsilon(episode)

			states_per_episode.append(episode_state)
			rewards_per_episodes.append(np.sum(total_reward))
			print(f"Episode {episode + 1}: Total Reward: {np.sum(total_reward)}")

			# decay the exploration rate proportional to its current value

			end_time = timeit.default_timer()

		for i, data in enumerate(rewards_per_episodes):
			print(data)
			logging.info(episode_state[i])

		print("-------------------------------------------------------")
		print(f"{data[0].name} test done in {end_time - start_time}")

if __name__ == '__main__':
	main()
