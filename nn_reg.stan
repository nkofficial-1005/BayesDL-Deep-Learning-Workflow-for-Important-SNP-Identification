functions {
  vector nn_predict(matrix x, matrix d_t_h, matrix[] h_t_h, vector h_t_d, row_vector[] hidden_bias, real y_bias) {
    int N = rows(x);
    int n_H = cols(d_t_h);
    int H = size(hidden_bias);
    matrix[N, n_H] hidden_layers[H];
    vector[N] output_layer;
    vector[N] ones = rep_vector(1., N);

    hidden_layers[1] = tanh(x * d_t_h + ones * hidden_bias[1]);
    for(h in 2:H) {
      hidden_layers[h] = tanh(hidden_layers[h-1] * h_t_h[h - 1] + ones * hidden_bias[h]);
    }
    output_layer = hidden_layers[H] * h_t_d + y_bias;
    return(output_layer);
  }		
}

data {
  int N; // Number of training samples
  int P; // Number of predictors (features)
  matrix[N, P] x; // Feature data
  vector[N] y; // Outcome
  int H; // Number of hidden layers
  int n_H; // Number of nodes per layer (All get the same)

  int N_test; // Number of test samples
  matrix[N_test, P] x_test; // Test predictors
}

transformed data {
}

parameters {
  matrix[P, n_H] data_to_hidden_weights; // Data -> Hidden 1
  matrix[n_H, n_H] hidden_to_hidden_weights[H - 1]; // Hidden[t] -> Hidden[t+1]
  vector[n_H] hidden_to_data_weights;
  // ordered[n_H] hidden_bias[H]; // Use ordered if using NUTS
  row_vector[n_H] hidden_bias[H]; // Hidden layer biases
  real y_bias; // Bias. 
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] output_layer;

  output_layer = nn_predict(x,
                            data_to_hidden_weights,
                            hidden_to_hidden_weights,
                            hidden_to_data_weights,
                            hidden_bias,
                            y_bias);

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
  y_bias ~ std_normal();

  sigma ~ std_normal();


  y ~ normal(output_layer, sigma);
}

generated quantities {
  vector[N_test] output_test = nn_predict(x_test,
				          data_to_hidden_weights,
				          hidden_to_hidden_weights,
				          hidden_to_data_weights,
				          hidden_bias,
				          y_bias);
  vector[N_test] output_test_rng;
  for(n in 1:N_test) {
    output_test_rng[n] = normal_rng(output_test[n], sigma);
  }
}