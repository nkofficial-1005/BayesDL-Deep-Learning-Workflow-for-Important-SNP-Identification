functions {
  vector[] nn_predict(matrix x, matrix d_t_h, matrix[] h_t_h, matrix h_t_d, row_vector[] hidden_bias, row_vector y_bias) {
    int N = rows(x);
    int n_H = cols(d_t_h);
    int H = size(hidden_bias);
    int num_labels = cols(y_bias) + 1;
    matrix[N, n_H] hidden_layers[H];
    vector[num_labels] output_layer_logit[N];
    vector[N] ones = rep_vector(1., N);

    hidden_layers[1] = tanh(x * d_t_h + ones * hidden_bias[1]);
    
    for(h in 2:H) {
      hidden_layers[h] = tanh(hidden_layers[h-1] * h_t_h[h - 1] + ones * hidden_bias[h]);
    }
    for(n in 1:N) {
      output_layer_logit[n, 1] = 0.0;
      output_layer_logit[n, 2:num_labels] = (hidden_layers[H, n] * h_t_d + y_bias)';
    }
    return(output_layer_logit);
  }
}

data {
  int N; // Number of training samples
  int P; // Number of predictors (features)
  matrix[N, P] x; // Feature data
  int labels[N]; // Outcome labels
  int H; // Number of hidden layers
  int n_H; // Number of nodes per layer (All get the same)

  int N_test; // Number of test samples
  matrix[N_test, P] x_test; // Test predictors
}

transformed data {
  int num_labels = max(labels); // How many labels are there
}

parameters {
  matrix[P, n_H] data_to_hidden_weights; // Data -> Hidden 1
  matrix[n_H, n_H] hidden_to_hidden_weights[H - 1]; // Hidden[t] -> Hidden[t+1]
  matrix[n_H, num_labels - 1] hidden_to_data_weights; // Hidden[T] -> Labels. Base class gets 0.
  // ordered[n_H] hidden_bias[H]; // Use ordered if using NUTS
  row_vector[n_H] hidden_bias[H]; // Hidden layer biases
  row_vector[num_labels - 1] labels_bias; // Labels biases. Base class gets 0.
}

transformed parameters {
  vector[num_labels] output_layer_logit[N]; // Predicted output layer logits

  output_layer_logit = nn_predict(x,
                                  data_to_hidden_weights,
                                  hidden_to_hidden_weights,
                                  hidden_to_data_weights,
                                  hidden_bias,
                                  labels_bias);

}

model {
  // Priors
  to_vector(data_to_hidden_weights) ~ std_normal();

  for(h in 1:(H-1)) {
    to_vector(hidden_to_hidden_weights[h]) ~ std_normal();
  }

  to_vector(hidden_to_data_weights) ~ std_normal();

  for(h in 1:H) {
    to_vector(hidden_bias[h]) ~ std_normal();
  }
  labels_bias ~ std_normal();

  for(n in 1:N) { // Likelihood
    labels[n] ~ categorical_logit(output_layer_logit[n]);
  }
}

generated quantities {
  vector[num_labels] output_layer_logit_test[N_test] = nn_predict(x_test,
							   data_to_hidden_weights,
							   hidden_to_hidden_weights,
							   hidden_to_data_weights,
							   hidden_bias,
							   labels_bias);
  matrix[N_test, num_labels] output_test;
  for(n in 1:N_test) {
    output_test[n] = softmax(output_layer_logit_test[n])';
  }
  
}
