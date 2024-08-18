import numpy as np

class QLearningAgent:
	def __init__(self, action_space, state_space, learning_rate, discount_factor, epsilon, epsilon_decay, epsilon_min, epsilon_max, precision):
		self.action_space = action_space
		self.state_space = state_space
		self.learning_rate = learning_rate
		self.discount_factor = discount_factor

		self.epsilon_rate = epsilon # exploration_rate
		self.epsilon_decay = epsilon_decay
		self.epsilon_min = epsilon_min
		self.epsilon_max = epsilon_max

		_total_discrete_actions = np.prod(self.action_space.nvec)
		_state_space_size = self.state_space.shape[0]
		self.q_table = np.zeros((_state_space_size, _total_discrete_actions), dtype=precision)

	def choose_action(self, state):
		if np.random.rand() < self.epsilon_rate:
			sampled_action = self.action_space.sample() # Random actions for all parameters and base stations
			action = sampled_action - 1

		else:
			action = np.argmax(self.q_table[state, :]) # best-known action using argmax

		return action


	def learn(self, state, action, reward, next_state):
		best_next_action = np.argmax(self.q_table[next_state]) # best-known action for next state
		td_target = reward + self.discount_factor * self.q_table[next_state][best_next_action]

		td_error = td_target - self.q_table[state][action]
		self.q_table[state][action] += self.learning_rate * td_error

		# from deeplizard Tutorial (this is awesome btw)
		# --------------------------------------------------
		# new_q_value = weighted-sum of old value and new value
		# q-new = (1 - alpha)(q-old) + alpha(reward + discount_factor * agmax(actions at q-new-state))

	def update_epsilon(self, episode_number):
		#self.epsilon *= self.epsilon_decay
		self.epsilon_rate = self.epsilon_min + (self.epsilon_max - self.epsilon_min) * np.exp(-self.epsilon_decay * episode_number)

	@property
	def status(self):
		return { 'epsilon' : self.epsilon_rate, 'q_table' : self.q_table }