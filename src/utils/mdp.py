import time, timeit
import subprocess
import numpy as np
from utils.common import *
from gym import spaces
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

		env_capture  = {}
		for tx_idx in range(self.num_transmitters):
			env_capture[tx_idx] = {}

		start_time = timeit.default_timer()

		if Path(self.model_bin).is_file():
			print(f"Model path ok: {self.model_bin}")

		print(f"Command\n{[self.model_bin] + simulation_parameters}")

		process = subprocess.Popen([self.model_bin] + simulation_parameters, \
			stdout=subprocess.PIPE, \
			stderr=subprocess.PIPE, bufsize=1)

		while True:
			stdout_line_str = process.stdout.readline().decode('utf-8')

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
			# e.g. "timeslot: 0  power: 2.00 alpha: 0.00 placement: 750,200  cow_id: 2 sta_id: 5  sinr: 13.95"
			elif not stdout_line_str.startswith("timeslot: "):
				self.logger.info(f"{stdout_line_str.rstrip()}")
				continue
			# end if

			try:
				stdout_line_str = stdout_line_str.rstrip()

				if stdout_line_str:
					self.logger.info(f"subprocess: {stdout_line_str}")

					tokens = stdout_line_str.split()

					tx_idx = None
					sta_id = None
					sinr = None

					for i in range(len(tokens)):
						if tokens[i] == 'cow_id:':
							tx_idx = int(tokens[i + 1])
						elif tokens[i] == 'sta_id:':
							sta_id = int(tokens[i + 1])
						elif tokens[i] == 'sinr:':
							sinr = float(tokens[i + 1])

					# Assign sinr to the correct cow_id and sta_id
					if sta_id not in env_capture[tx_idx]:
						env_capture[tx_idx][sta_id] = sinr

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

		self.logger.info(f"get_output() {env_capture}")
		return env_capture
	# end def

	def prime_subprocess(self):

		_location = {
			'tx_station' : np.transpose([self.arguments['transmitter_location_x'], self.arguments['transmitter_location_y']]),
			'rx_station' : self.receiver_locations
		}


		# validate
		sim_args = [
			'--timeslot',				str(self.timeslot),\
			'--frequency',				str(self.frequency),\
			'--bandwidth',				str(self.bandwidth),\
			'--symbolrate',				str(self.symbolrate),\
			'--blockpersymbol',			str(self.blockpersymbol),\
			'--height',					str(self.antenna_height),\
			'--ms_grx',					str(self.ms_grx),\
			'--system_noise',			str(self.system_noise),\
			'--base_stations',			str(self.num_transmitters),\
			'--mobile_stations',		str(self.num_receivers),\
			#'--timeslots',			    '1',\
			'--slimit',			        str(self.receiver_snr_limit),\
			'--antenna_dims',                      ','.join([str(s) for s in self.antenna_size]),\
			'--base_station_power_range_dBm',      ','.join([str(s) for s in self.limit_power]),\
			'--base_station_scan_angle_range_deg', ','.join([str(s) for s in self.limit_scana]),\
			'--antenna_txpower' ,                  ','.join([str(s) for s in self.arguments['antenna_txpower']]),\
			'--scan_angle'      ,                  ','.join([str(s) for s in self.arguments['antenna_scan_angle']]),\
			'--antenna_spacing',                   ','.join([str(s) for s in self.arguments['antenna_spacing']]),\
			'--base_station_antenna_counts',       ','.join([str(s) for s in self.arguments['antenna_counts']]),\
			#'--base_station_antenna_counts',       ','.join([str(s) for s in self.antenna_panels]),\
			'--base_station_theta_c',              ','.join([str(s) for s in self.arguments['transmitter_direction']]),\
			'--base_station_location',             ' '.join([f"{x},{y}" for x, y in _location['tx_station']]),\
			'--mobile_station_location',           ' '.join([f"{x},{y}" for x, y in _location['rx_station']]),\
			'--ms_selection',                      ','.join([str(s) for s in self.arguments['receiver_selection_index']])
		]

		return sim_args
	# end def

	def update_args(self, actions_list):
		index = 0

		for key, parameter_set in self.arguments:
			param_info = self.arguments_metadata[key]
			param_type = param_info['type']
			param_range = param_info['range']

			range_span = (param_range[1] - param_range[0]) if param_type == 'continuous' else 0

			# tx is transmitter
			for tx_idx in range(self.num_transmitters):
				discrete_action = actions_list[index]
				action_per_parameter_per_tx = discrete_action - 1  # -1, 0, +1

				if param_type == 'continuous':
					# Calculate delta based on the range and action

					#randomness_factor = np.random.uniform(param_range[0], param_range[1])
					delta = action_per_parameter_per_tx * (range_span / 100) # * randomness_factor

				else:# 'discrete'
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

		self.update_args(action_list)

		args = self.prime_subprocess()
		self.last_state = self.get_output(args)

		sinr_list = []
		for tx_idx in range(self.num_transmitters):
			sinr_list.append(self.last_state[tx_idx][self.arguments['receiver_selection_index'][tx_idx]])

		reward, done = self.compute_reward(sinr_list)
		return self.get_state(), reward, done
	# end def

	# provides the Agent with the current state of the environment and make next decisions
	def get_state(self):
		state = np.array([[tx_id for tx_id in self.arguments[key]] for key in self.arguments])
		return state
	# end def

	# capture the state of the system from the last call
	@property
	def capture(self):
		return self.last_state

	def freeplay(self):
		args = self.prime_subprocess()

		self.last_state = self.get_output(args)

		print(self.last_state)

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
	def reset(self, timeslot = None, json_data = None):
		# either zero out the simulation or start with JSON init params (for timeslot #0, #1, #2, #3 ..)
		self.arguments.clear()
		self.timeslot = timeslot

		if not json_data:
			self.arguments.update({ # this can be randomized rather than use the json file ... between episodes
				'antenna_txpower'          : [0.0] * self.num_transmitters,
				'antenna_scan_angle'       : [0.0] * self.num_transmitters,
				'antenna_spacing'          : [0]   * self.num_transmitters,
				'antenna_counts'           : [1]   * self.num_transmitters,
				'transmitter_direction'    : [0.0] * self.num_transmitters,
				'transmitter_location_x'   : [0.0] * self.num_transmitters,
				'transmitter_location_y'   : [0.0] * self.num_transmitters,
				'receiver_selection_index' : [0]   * self.num_transmitters
			})

		else:
			_transmitter_locations = np.transpose(json_data['base_station_location'])
			self.arguments.update({
				'antenna_txpower'          : [float(txpower) for txpower in json_data['base_station_power_dBm_lut'][timeslot]],
				'antenna_scan_angle'       : [float(angle) for angle in json_data['base_station_scan_alpha_deg_lut'][timeslot]],
				'antenna_spacing'          : [float(spacing) for spacing in json_data['antenna_spacing']],
				'antenna_counts'           : [int(f) for f in json_data['base_station_antenna_counts']],
				'transmitter_direction'    : [float(f) for f in json_data['base_station_theta_c']],
				'transmitter_location_x'   : [float(coordinates) for coordinates in _transmitter_locations[0]],
				'transmitter_location_y'   : [float(coordinates) for coordinates in _transmitter_locations[1]],
				'receiver_selection_index' : [int(node_idx) - 1 for node_idx in json_data['base_to_mobile_station_id_selection_lut'][timeslot]],
			})

		return self.get_state()
	# end def

	# for advanced timeslots
	def get_mask(self):
		mask = []
		for value in self.arguments_metadata:
			if value['const_in_advance']:
				temp_arr = [0] * self.num_transmitters
			else:
				temp_arr = [1] * self.num_transmitters
			mask.extend(temp_arr)

		return mask

	@property
	def num_parameters(self):
		return len(self.arguments_metadata) * self.num_transmitters

	def __init__(self, binary, limit_sinr, debug, log_handle, json_data = None):

		global global_c_speed

		# validate the executatable
		self.model_bin         = binary
		self.limit_reward      = limit_sinr
		self.logger            = log_handle
		self.debug_flag        = debug
		self.timeout           = None
		self.last_state        = None
		self.timeslot          = None
		self.frequency         = float(json_data['frequency'])
		self.bandwidth         = 20e6 if not json_data else float(json_data['bandwidth'])
		self.symbolrate        = int(json_data['symbolrate'])
		self.blockpersymbol    = int(json_data['blockpersymbol'])
		self.antenna_height    = int(json_data['antenna_height'])
		self.antenna_size      = [float(s) for s in json_data['antenna_dims']]
		#self.antenna_panels    = [float(s) for s in json_data['base_station_antenna_counts']]
		self.ms_grx            = int(json_data['ms_grx_dB'])
		self.system_noise      = 5  if not json_data else float(json_data['system_noise_dB'])
		self.num_transmitters  = 3  if not json_data else int(json_data['base_stations'])
		self.num_receivers     = 15 if not json_data else int(json_data['mobile_stations'])

		self.receiver_locations = json_data['mobile_station_location']
		self.receiver_snr_limit = 24 if not json_data else float(json_data['sinr_limit_dB'])

		self.limit_power             = [-30.0, 30.0] if not json_data else [float(f) for f in json_data['base_station_power_range_dBm']]
		self.limit_scana             = [-90.0, 90.0] if not json_data else [float(f) for f in json_data['base_station_scan_angle_range_deg']]
		self.limit_antenna_count     = [1, 6] if not json_data else [float(f) for f in json_data['antenna_count_limit']]
		self.limit_antenna_direction = [0, 360] if not json_data else [float(f) for f in json_data['base_station_theta_c_lim_deg']]

		# discrete places where the clients are placed
		self.limit_tx_location_range = [[ 0, 0 ], [ 1000, 200 ]] if not json_data else \
			[[float(c) for c in coordinates] for coordinates in json_data['base_station_field_range_meters']]


		#self.limit_rx_location_range = [[x, y] for y in self.limit_rx_location_range[1] for x in self.limit_rx_location_range[0]]

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
		_arguments_bounds = [[], []]

		for key, value in self.arguments_metadata.items():
			_range = value['range']
			_arguments_bounds[0].extend([_range[0]] * self.num_transmitters)
			_arguments_bounds[1].extend([_range[1]] * self.num_transmitters)

		# Determine the action space. For each parameter, we'll have 3 possible actions: -1 (decrease), 0 (no change), +1 (increase)
		#self.action_space = spaces.Box(low=np.array([-1] * _num_parameters), high=np.array([1] * _num_parameters), dtype=np.float32)
		self.observation_space = spaces.Box(low=np.array(_arguments_bounds[0]), high=np.array(_arguments_bounds[1]), dtype=np.float32)
	# end def
