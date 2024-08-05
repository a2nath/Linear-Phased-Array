import time, timeit
import subprocess
import numpy as np
from utils.common import *
import gym
from collections import OrderedDict

import pdb

class MDPEnv:

	def getrng_int(self, key, range=[0, 1], size=(1,1)):

		#if np.random.rand() < self.epsilon:
		#	# Exploration: Choose a random value in the range
		#	rand = np.random.randint(range[0], range[1], size)
		#else:
		#	# Exploitation: Use the known best value with some randomness
		#	rand = self.arguments[key] + np.random.randint(-1, 2, size)


		rand = np.random.randint(self.arguments_metadata[key][0], self.arguments_metadata[key][0])

		new_value = np.clip(rand, range[0], range[1])

		return new_value
	# end def

	def getrng_float(self, key, range=[0.0, 1.0], samples=1):

		if np.random.rand() < self.epsilon:
			# Exploration: Choose a random value in the range
			rand = np.random.uniform(range[0], range[1], samples)
		else:
			# Exploitation: Use the known best value with some randomness
			rand = self.arguments[key] + np.random.uniform(-0.1, 0.1, samples)

		new_value = np.clip(rand, range[0], range[1])
		return new_value
	# end def

	#def get_next_intval(self, key, range=[0, 1]):
	#
	#

	def get_next_state(self):
		#states = get_states(self.arguments)
		#Q = np.zeros((1, states))
		return

	def get_output(self, simulation_parameters):
		global global_timeout
		sinr_data = []
		start_time = timeit.default_timer()

		if Path(self.model_bin).is_file():
			print("model path ok")

		print([self.model_bin] + self.model)

		process = subprocess.Popen([self.model_bin] + simulation_parameters, \
			stdout=subprocess.PIPE, \
			stderr=subprocess.PIPE, bufsize=1)

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
				# end while
			# end if

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
			#end try

			if global_timeout is not None and float(timeit.default_timer() - start_time) > float(global_timeout):
				process.kill()
				outs, errs = process.communicate(timeout=global_timeout)
				print("----------------------KILLED---------------------------")
				print(f"Process timed out after {global_timeout} seconds\nouts = {outs.decode('utf-8').rstrip()}\nerrs = {errs.decode('utf-8').rstrip()}")
		# end while

		self.logger.info(sinr_data)

		print(f"sinr_data {sinr_data}")
		return sinr_data
	# end def

	def prime_subprocess(self):

		_location = {
			'tx_station' : np.transpose([self.arguments['transmitter_location_x'], self.arguments['transmitter_location_y']]),
			'rx_station' : np.transpose(self.limit_rx_location_range)[self.arguments['receiver_selection_index']]
		}


		# validate
		sim_args = [
			'--frequency',				str(self.frequency),\
			'--bandwidth',				str(self.bandwidth),\
			'--symbolrate',				str(self.symbolrate),\
			'--blockpersymbol',			str(self.blockpersymbol),\
			'--height',					str(self.antenna_height),\
			'--ms_grx',					str(self.ms_grx),\
			'--system_noise',			str(self.system_noise),\
			'--mobile_stations',		str(self.num_transmitters),\
			'--base_stations',			str(self.num_receivers),\
			'--timeslots',			    '1',\
			'--slimit',			        str(self.limit_sinr),\
			'--antenna_dims',                      ','.join([str(s) for s in self.antenna_size]),\
			'--base_station_power_range_dBm',      ','.join([str(s) for s in self.limit_power]),\
			'--base_station_scan_angle_range_deg', ','.join([str(s) for s in self.limit_scana]),\
			'--antenna_txpower' ,                  ','.join([str(s) for s in self.arguments['antenna_txpower']]),\
			'--scan_angle'      ,                  ','.join([str(s) for s in self.arguments['antenna_scan_angle']]),\
			'--antenna_spacing',                   ','.join([str(s) for s in self.arguments['antenna_spacing']]),\
			'--base_station_antenna_counts',       ','.join([str(s) for s in self.arguments['antenna_counts']]),\
			'--base_station_theta_c',              ','.join([str(s) for s in self.arguments['transmitter_direction'].value]),\
			'--base_station_location',             ','.join([str(s) for s in _location['tx_station']]),\
			'--mobile_station_location',           ','.join([str(s) for s in _location['rx_station']]),\
			'--ms_selection',                      ','.join([str(s) for s in self.arguments['receiver_selection_index']])
		]

		return sim_args
	# end def

	def update_args(self, actions_list):
		index = 0;

		for key, parameter_set in self.arguments:
			param_info = self.arguments_metadata[key]
			param_type = param_info['type']
			param_range = param_info['range']

			range_span = (param_range[1] - param_range[0]) if param_type == 'continuous' else 0

			# tx is transmitter
			for tx_idx in range(self.num_transmitters):
				action_per_parameter_per_tx = actions_list[index] - 1  # -1, 0, +1

				if param_type == 'continuous':
					# Calculate delta based on the range and action

					#randomness_factor = np.random.uniform(param_range[0], param_range[1])
					delta = action_per_parameter_per_tx * (range_span / 100) # * randomness_factor

				else:# 'discrete'
					action_per_parameter_per_tx = round(action_per_parameter_per_tx)
					# For discrete, ensure the action results in a valid value
					delta = action_per_parameter_per_tx
					#delta = action_per_parameter_per_tx * (range_span / 100) # * randomness_factor
					#new_index = self.action_index[key] + action_per_parameter_per_tx
					#new_index = np.clip(new_index, param_range)
					#new_value = param_range[new_index]

				# Clip the new value within the parameter's range
				new_value = np.clip(parameter_set[tx_idx] + delta, param_range[0], param_range[1])
				parameter_set[tx_idx] = new_value

				index += 1
	# end def

	def step(self, action_list):

		self.update_args(action_list);

		args = self.prime_subprocess()
		sinr_values = self.get_output(args)
		reward, done = self.compute_reward(sinr_values)
		return self.get_state(), reward, done
	# end def

	# provides the Agent with the current state of the environment and make next decisions
	def get_state(self):
		state = np.array([[tx_id for tx_id in self.arguments[key]] for key in self.arguments])
		return state
	# end def

	# Check if the episode is done based on the sinr_values
	def check_done(self, sinr_values):
		return sum(float(snr) >= self.limit_reward for snr in sinr_values)
	# end def

	# Compute reward based on the sinr_values
	def compute_reward(self, sinr_values):
		_num_signal_requirement_met = self.check_done(sinr_values);

		if _num_signal_requirement_met >= self.num_transmitters:
			return self.num_transmitters, True

		elif _num_signal_requirement_met < len(sinr_values):
			return _num_signal_requirement_met * 0.5, False # revisit this if needed

		elif _num_signal_requirement_met == 0:
			return -10, False
	# end def

	@property
	def results(self):
		return self.arguments

	# timslot_num > 0, restrict the parameters
	def advance_simulation(self):
		for index, (key, param_info) in enumerate(self.arguments_metadata.items()):
			param_range = param_info['range']

			for tx_idx in range(self.num_transmitters):
				# in the next timeslot p_tx can be +/- 1 dB
				if key == 'antenna_txpower':
					current_power = self.arguments['antenna_txpower'][tx_idx]
					self.observation_space.low [index + tx_idx] = max(current_power - 1.0, param_range[0])
					self.observation_space.high[index + tx_idx] = min(current_power + 1.0, param_range[1])

				# keep the parameters const across timeslots once found the best
				elif param_info['const_in_advance']:
					current_value = self.arguments[key][tx_idx]
					self.observation_space.low [index + tx_idx] = current_value
					self.observation_space.high[index + tx_idx] = current_value
	# end def

	# use between episodes, <<< not timeslots! >>>
	def reset(self, timeslot, json_data):
		# either zero out the simulation or start with JSON init params (for timeslot #0, #1, #2, #3 ..)
		self.arguments.clear()

		if not json_data:

			self.arguments.update({
				'antenna_txpower'          : [0.0]*self.num_transmitters,
				'antenna_scan_angle'       : [0.0]*self.num_transmitters,
				'antenna_spacing'          : [0]  *self.num_transmitters,
				'antenna_counts'           : [1]  *self.num_transmitters,
				'transmitter_direction'    : [0.0]*self.num_transmitters,
				'transmitter_location_x'   : [0.0]*self.num_transmitters,
				'transmitter_location_y'   : [0.0]*self.num_transmitters,
				'receiver_selection_index' : [0]  *self.num_transmitters
			})

		else:
			j_antenna_txpower = json_data['base_station_power_dBm_lut'][timeslot]
			j_scan_angle      = json_data['base_station_scan_alpha_deg_lut'][timeslot]
			j_ms_selection    = json_data['base_to_mobile_station_id_selection_lut'][timeslot]

			_transmitter_locations = np.transpose(json_data['base_station_location'])

			self.arguments.update({
				'antenna_txpower'          : [float(txpower) for txpower in j_antenna_txpower],
				'antenna_scan_angle'       : [float(angle) for angle in j_scan_angle],
				'antenna_spacing'          : [float(spacing) for spacing in json_data['antenna_spacing']],
				'antenna_counts'           : [int(f) for f in json_data['base_station_antenna_counts']],
				'transmitter_direction'    : [float(f) for f in json_data['base_station_theta_c']],
				'transmitter_location_x'   : [float(coordinates) for coordinates in _transmitter_locations[0]],
				'transmitter_location_y'   : [float(coordinates) for coordinates in _transmitter_locations[1]],
				'receiver_selection_index' : [int(node_idx) for node_idx in j_ms_selection],
			})

		return self.get_state()
	# end def

	def __init__(self, binary, limit_sinr, debug, log_handle, json_data = None):

		global global_c_speed

		# validate the executatable
		self.model_bin         = binary
		self.limit_reward      = limit_sinr
		self.logger            = log_handle
		self.debug_flag        = debug
		self.timeout           = None
		self.frequency         = float(json_data['frequency']),
		self.bandwidth         = 20e6 if not json_data else float(json_data['bandwidth']),
		self.symbolrate        = int(json_data['symbolrate']),
		self.blockpersymbol    = int(json_data['blockpersymbol']),
		self.antenna_height    = int(json_data['antenna_height']),
		self.antenna_size      = int(json_data['antenna_dims']),
		self.ms_grx            = int(json_data['ms_grx_dB']),
		self.system_noise      = 5  if not json_data else float(json_data['system_noise_dB']),
		self.num_transmitters  = 3  if not json_data else int(json_data['mobile_stations']),
		self.num_receivers     = 15 if not json_data else int(json_data['base_stations']),
		self.limit_sinr        = 24 if not json_data else float(json_data['sinr_limit_dB']),

		self.limit_power             = [-30.0, 30.0] if not json_data else [float(f) for f in json_data['base_station_power_range_dBm']]
		self.limit_scana             = [-90.0, 90.0] if not json_data else [float(f) for f in json_data['base_station_scan_angle_range_deg']]
		self.limit_antenna_count     = [1, 6] if not json_data else [float(f) for f in json_data['antenna_count_limit']]
		self.limit_antenna_direction = [0, 360] if not json_data else [float(f) for f in json_data['base_station_theta_c_lim_deg']]

		# discrete places where the clients are placed
		self.limit_tx_location_range = [[ 0, 0 ], [ 1000, 200 ]] if not json_data else [float(f) for f in json_data['base_station_field_range_meters']]
		self.limit_rx_location_range = [list(range(200, 800, 50)), list(range(450, 650, 50))] if not json_data else [float(f) for f in json_data['mobile_station_field_range_meters']]

		_limit_tx_location_range = np.transpose(self.limit_tx_location_range)
		#self.receiver_list = list(range(self.num_receivers))

		#self.action_idx = {}
		self.arguments = OrderedDict()

		# some parameters are discrete, or constant across timeslots
		self.arguments_metadata = OrderedDict()
		self.arguments_metadata.update({
			'antenna_txpower' : {'type': 'continuous', 'const_in_advance' : True, 'range': self.limit_power},
			'antenna_scan_angle' : {'type': 'continuous', 'const_in_advance' : False, 'range': self.limit_scana},
			'antenna_spacing' : {'type': 'continuous', 'const_in_advance' : True, 'range': [0, 2 * global_c_speed/self.frequency]},
			'antenna_counts' : {'type': 'discrete', 'const_in_advance' : False, 'range': self.limit_antenna_count},
			'transmitter_direction' : {'type': 'continuous', 'const_in_advance' : True, 'range': self.limit_antenna_direction},
			'transmitter_location_x' : {'type': 'continuous', 'const_in_advance' : True, 'range':_limit_tx_location_range[0]},
			'transmitter_location_y' : {'type': 'continuous', 'const_in_advance' : True, 'range':_limit_tx_location_range[1]},
			'receiver_selection_index' : {'type': 'discrete', 'const_in_advance' : False, 'range': [0, self.num_receivers - 1]}
		})

		# Automatically scale the number of actions
		_num_parameters = len(self.arguments_metadata) * self.num_transmitters
		_arguments_bounds = [[], []]

		for key, value in self.arguments_metadata.items():
			_range = value['range']
			_arguments_bounds[0].extend([_range[0]] * self.num_transmitters)
			_arguments_bounds[1].extend([_range[1]] * self.num_transmitters)

		# Determine the action space. For each parameter, we'll have 3 possible actions: -1 (decrease), 0 (no change), +1 (increase)
		self.action_space = spaces.Box(low=np.array([-1] * _num_parameters), high=np.array([1] * _num_parameters), dtype=np.float32)
		self.observation_space = spaces.Box(low=np.array(_arguments_bounds[0]), high=np.array(_arguments_bounds[1]), dtype=np.float32)
	# end def

class QLearningAgent:
	def __init__(self, action_space, state_space, learning_rate=0.1, discount_factor=0.99, epsilon=0.1, epsilon_decay=0.99):
		self.action_space = action_space
		self.state_space = state_space
		self.learning_rate = learning_rate
		self.discount_factor = discount_factor
		self.epsilon = epsilon
		self.epsilon_decay = epsilon_decay
		self.q_table = np.zeros((state_space, action_space))

	def choose_action(self, state):

		if np.random.rand() < self.epsilon:
			action = self.action_space.sample()  # Random actions for all parameters and base stations
		else:
			state_index = self.state_to_index(state)
			action = np.argmax(self.q_table[state_index]) # best-known action using argmax

		return action


	def learn(self, state, action, reward, next_state):
		best_next_action = np.argmax(self.q_table[next_state]) # best-known action for next state
		td_target = reward + self.discount_factor * self.q_table[next_state][best_next_action]
		td_error = td_target - self.q_table[state][action]
		self.q_table[state][action] += self.learning_rate * td_error

	def update_epsilon(self):
		self.epsilon *= self.epsilon_decay
