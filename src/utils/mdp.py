from argument import *
from common import *
import time, timeit
import subprocess
import numpy as np

class Argument:
	# val        : value of the argument
	# mutable    : can be changed between steps RL
	# parameters : counts as 1 or more parameter(s)
	def __init__(self, val, mutable = False, parameters = 1):
		self.val = val
		self.mutable = mutable
		self.num_parameters = parameters


class MDP:

	def genrng_int(self, key, range=[0, 1], size=(1,1)):

		if np.random.rand() < self.epsilon:
		    # Exploration: Choose a random value in the range
		    rand = np.random.randint(range[0], range[1], size)
		else:
		    # Exploitation: Use the known best value with some randomness
		    rand = self.sim_args[key] + np.random.randint(-1, 2, size)

		new_value = np.clip(rand, range[0], range[1])

		return new_value


	def getrng_float(self, key, range=[0.0, 1.0], samples=1):

		if np.random.rand() < self.epsilon:
		    # Exploration: Choose a random value in the range
		    rand = np.random.uniform(range[0], range[1], samples)
		else:
		    # Exploitation: Use the known best value with some randomness
		    rand = self.sim_args[key] + np.random.uniform(-0.1, 0.1, samples)

		new_value = np.clip(rand, range[0], range[1])
		return new_value

	#def get_next_intval(self, key, range=[0, 1]):
	#
	#

	def get_next_state(self):
		states = get_states(self.arguments)
		Q = np.zeros((1, states))


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

	def update_args(self, key):
		sle

	def init_args(self):

		# validate
		j_antenna_txpower = self.arguments['base_station_power_dBm_lut']
		j_scan_angle      = self.arguments['base_station_scan_alpha_deg_lut']
		j_ms_selection    = self.arguments['base_to_mobile_station_id_selection_lut']

		timeslots = 1

		if len(j_antenna_txpower) == len(j_scan_angle) and len(j_scan_angle) == len(j_ms_selection):
			timeslots = len(j_antenna_txpowre)
		else:
			raise ValueError("The number of timeslots is mismatch between power, scan, ms:bs id mapping.")

		self.sim_args = {
			'--frequency',				str(self.arguments['frequency']),\
			'--bandwidth',				str(self.arguments['bandwidth']),\
			'--symbolrate',				str(self.arguments['symbolrate']),\
			'--blockpersymbol',			str(self.arguments['blockpersymbol']),\
			'--height',					str(self.arguments['height']),\
			'--ms_grx',					str(self.arguments['ms_grx']),\
			'--system_noise',			str(self.arguments['system_noise']),\
			'--mobile_stations',		str(self.arguments['mobile_stations']),\
			'--base_stations',			str(self.arguments['base_stations']),\
			'--timeslots',			    str(self.arguments['timeslots']),\
			'--slimit',			        str(self.arguments['slimit']),\

			# 3 parameters
			'--base_station_theta_c',	','.join(map(str, self.arguments['base_station_theta_c'].val)),\
			# 6 parameters
			'--base_station_location',		convert_lut(self.arguments['base_station_location']),\

			'--mobile_station_location',	convert_lut(self.arguments['mobile_station_location']),\
			# 3 parameters
			'--base_station_antenna_counts',			','.join(map(str, self.arguments['base_station_antenna_counts'])),\

			'--base_station_power_range_dBm',			','.join(map(str, self.arguments['base_station_power_range_dBm'])),\
			'--base_station_scan_angle_range_deg',		','.join(map(str, self.arguments['base_station_scan_angle_range_deg'])),\

			# 3 parameters
			'--antenna_spacing',		','.join(map(str, self.arguments['antenna_spacing'])),\

			'--antenna_dims',			','.join(map(str, self.arguments['antenna_dims'])),\

			# 3 parameters
			'--antenna_txpower',            convert_lut(j_antenna_txpower),\

			# 3 parameters
			'--scan_angle',	                convert_lut(j_scan_angle),\

			# 3 parameters
			'--ms_selection',	            reparse_binding(j_ms_selection['base_to_mobile_station_id_selection_lut']),\
		}


		for key, argument in self.arguments.items():
			if argument.mutable == True:
				self.num_parameters += argument.num_parameters


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

		#self.Q[] = average_sinr / len(output)

		#print(sirn_data)
		#last_action

	def build_args(self, arguments, jdata: list):

		# validate
		j_antenna_txpower = jdata['base_station_power_dBm_lut']
		j_scan_angle      = jdata['base_station_scan_alpha_deg_lut']
		j_ms_selection    = jdata['base_to_mobile_station_id_selection_lut']

		timeslots = 1

		if len(j_antenna_txpower) == len(j_scan_angle) and len(j_scan_angle) == len(j_ms_selection):
			timeslots = len(j_antenna_txpowre)
		else:
			raise ValueError("The number of timeslots is mismatch between power, scan, ms:bs id mapping.")

		self.arguments['--frequency']       = Argument(float(jdata['frequency']), mutable=False, parameters=1)
		self.arguments['--bandwidth']       =	Argument(float(jdata['bandwidth']), mutable=False, parameters=1)
		self.arguments['--symbolrate']      = Argument(int(jdata['symbolrate']), mutable=False, parameters=1)
		self.arguments['--blockpersymbol']  = Argument(int(jdata['blockpersymbol']), mutable=False, parameters=1)
		self.arguments['--height']          = Argument(int(jdata['height']), mutable=False, parameters=1)
		self.arguments['--ms_grx']          = Argument(int(jdata['ms_grx_dB']), mutable=False, parameters=1)
		self.arguments['--system_noise']    = Argument(float(jdata['system_noise']), mutable=False, parameters=1)
		self.arguments['--mobile_stations'] = Argument(int(jdata['mobile_stations']), mutable=False, parameters=1)
		self.arguments['--base_stations']   = Argument(int(jdata['base_stations']), mutable=False, parameters=1)
		self.arguments['--timeslots']       = Argument(int(timeslots), mutable=False, parameters=1) # N
		self.arguments['--slimit']          = Argument(float(jdata['sinr_limit']), mutable=False, parameters=1)

		base_stations = self.arguments['--base_stations']

		# 3 parameters
		self.arguments['--base_station_theta_c']              = Argument(map(float, jdata['base_station_theta_c']), mutable=True, parameters=base_stations)
		# 6 parameters
		self.arguments['--base_station_location']             =	\
			Argument([map(float, coordinates) for coordinates in jdata['base_station_location']], mutable=True, parameters=base_stations)

		self.arguments['--mobile_station_location']           = \
			Argument([map(float, coordinates) for coordinates in jdata['mobile_station_location']], mutable=False, parameters=1)
		# 3 parameters
		self.arguments['--base_station_antenna_counts']       = Argument(int(jdata['base_station_antenna_counts']), mutable=True, parameters=base_stations)

		self.arguments['--base_station_power_range_dBm']      = Argument(map(float, jdata['base_station_power_range_dBm']), mutable=False, parameters=1)
		self.arguments['--base_station_scan_angle_range_deg'] =	Argument(map(float, jdata['base_station_scan_angle_range_deg']), mutable=False, parameters=1)

		# 3 parameters
		self.arguments['--antenna_spacing']  = Argument(map(float, jdata['antenna_spacing']), mutable=True, parameters=base_stations)

		self.arguments['--antenna_dims']     = Argument(map(float, jdata['antenna_dims']), mutable=False, parameters=1)

		# 3N parameters
		self.arguments['--antenna_txpower']  = Argument([map(float, coordinates) for coordinates in j_antenna_txpower], mutable=True, parameters=base_stations*timeslots)

		# 3N parameters
		self.arguments['--scan_angle']       = Argument([map(float, coordinates) for coordinates in j_scan_angle], mutable=True, parameters=base_stations*timeslots)

		# 3N parameters
		self.arguments['--ms_selection']     = Argument([map(int, coordinates) for coordinates in j_ms_selection], mutable=True, parameters=base_stations*timeslots)



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
		self.epsilon          = args.epsilon
		self.timeout          = None            # for subprocess
		self.sim_args         = None            # for passing the args to bin
		self.num_parameters   = None            # for size of RNG
		self.arguments        = {}

		self.arguments['frequency'] = Argument


		#else:
		#	raise FileNotFoundError(f"Executable file has a problem or does not exist {default_model_bin}")


