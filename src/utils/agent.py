import numpy as np

class LinearQFunction:
    def __init__(self, num_features, num_actions, alpha=0.1, default_q_value=0.0):
        self.num_features = num_features
        self.num_actions = num_actions
        self.alpha = alpha
        self.weights = np.full((num_actions, num_features), default_q_value, dtype=np.uint8)

    def get_q_value(self, feature_vector, action):
        return np.dot(self.weights[action], feature_vector)

    def update(self, feature_vector, action, delta):
        self.weights[action] += self.alpha * delta * feature_vector

class FeatureExtractor:
    def __init__(self, state_space_size, action_space_size):
        self.state_space_size = state_space_size
        self.action_space_size = action_space_size

    def extract(self, state, action):
        feature_vector = np.zeros(self.state_space_size * self.action_space_size, dtype=np.float32)
        state_index = state  # Assuming state is already an index or can be mapped to one
        action_index = action
        feature_index = state_index * self.action_space_size + action_index
        feature_vector[feature_index] = 1.0
        return feature_vector

class QLearningAgent:
	def __init__(self, action_space, state_space_size, learning_rate, discount_factor, epsilon, epsilon_decay, epsilon_min, epsilon_max, precision):
		self.action_space = action_space
		self.learning_rate = learning_rate
		self.discount_factor = discount_factor

		self.epsilon_rate = epsilon # exploration_rate
		self.epsilon_decay = epsilon_decay
		self.epsilon_min = epsilon_min
		self.epsilon_max = epsilon_max

		# Initialize function approximator and feature extractor
		self.num_actions = np.prod(action_space.nvec)
		self.feature_extractor = FeatureExtractor(state_space_size, self.num_actions)
		self.q_function = LinearQFunction(num_features=state_space_size * self.num_actions, num_actions=self.num_actions, alpha=learning_rate)

	def choose_action(self, state):
		if np.random.rand() < self.epsilon_rate:
			sampled_action = self.action_space.sample() # Random actions for all parameters and base stations
			action = sampled_action - 1
		else:
			#action = np.argmax(self.q_table[state, :]) # best-known action using argmax
			q_values = [self.q_function.get_q_value(self.feature_extractor.extract(state, action), action) for action in range(self.num_actions)]
			action = np.argmax(q_values)

		return action



	def learn(self, state, action, reward, next_state):
		# Compute target
		next_q_values = [self.q_function.get_q_value(self.feature_extractor.extract(next_state, a), a) for a in range(self.num_actions)]
		best_next_action = np.argmax(next_q_values)
		td_target = reward + self.discount_factor * next_q_values[best_next_action]

		# Update the Q-function using the TD error
		current_q_value = self.q_function.get_q_value(self.feature_extractor.extract(state, action), action)
		td_error = td_target - current_q_value
		self.q_function.update(self.feature_extractor.extract(state, action), action, td_error)




		# from deeplizard Tutorial (this is awesome btw)
		# --------------------------------------------------
		# new_q_value = weighted-sum of old value and new value
		# q-new = (1 - alpha)(q-old) + alpha(reward + discount_factor * agmax(actions at q-new-state))

	def update_epsilon(self, episode_number):
		#self.epsilon *= self.epsilon_decay
		self.epsilon_rate = self.epsilon_min + (self.epsilon_max - self.epsilon_min) * np.exp(-self.epsilon_decay * episode_number)

	@property
	def status(self):
		return { 'epsilon' : self.epsilon_rate, 'q_table' : self.q_function.weights }