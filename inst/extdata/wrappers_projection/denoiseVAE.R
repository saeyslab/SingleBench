
wrapper.projection.denoiseVAE <- WrapTool(
  name = 'denoiseVAE',
  type = 'projection',
  r_packages = c('torch'),
  fun.build_model =
    function(
      input,
      latent_dim,
      hidden_dim = 20,
      n_epochs = 20,
      learning_rate = 0.001,
      batch_size = 512
    ) {
      
      vae <- torch::nn_module(
        initialize = function(input_dim, hidden_dim, latent_dim) {
          self$encoder <- torch::nn_sequential(
            torch::nn_linear(input_dim, hidden_dim),
            torch::nn_relu(),
            torch::nn_linear(hidden_dim, latent_dim * 2) # mu, log_var
          )
          self$decoder <- torch::nn_sequential(
            torch::nn_linear(latent_dim, hidden_dim),
            torch::nn_relu(),
            torch::nn_linear(hidden_dim, input_dim),
            torch::nn_sigmoid()
          )
        },
        
        reparameterize = function(mu, log_var) {
          std <- torch::torch_exp(0.5 * log_var)
          eps <- torch::torch_randn_like(std)
          mu + eps * std
        },
        
        forward = function(x) {
          encoded <- self$encoder(x)
          mu <- encoded[, 1:(ncol(encoded) / 2)]
          log_var <- encoded[, (ncol(encoded) / 2 + 1):ncol(encoded)]
          z <- self$reparameterize(mu, log_var)
          decoded <- self$decoder(z)
          list(decoded, mu, log_var)
        }
      )
      
      input_dim <- ncol(input)
      
      model <- vae(input_dim, hidden_dim, latent_dim)
      
      loss_fn <- function(recon_x, x, mu, log_var) {
        recon_loss <- torch::nn_mse_loss()(recon_x, x)
        kld_loss <- -0.5 * torch::torch_sum(1 + log_var - mu^2 - torch::torch_exp(log_var))
        recon_loss + kld_loss
      }
      optimizer <- torch::optim_adam(model$parameters, lr = learning_rate)
      
      data_loader <- torch::dataloader(
        dataset = torch::tensor_dataset(torch::torch_tensor(input)),
        batch_size = batch_size,
        shuffle = TRUE
      )
      
      for (epoch in seq_len(n_epochs)) {
        total_loss <- 0
        for (batch in torch::enumerate(data_loader)) {
          batch_data <- batch[[1]]
          
          optimizer$zero_grad()
          outputs <- model(batch_data)
          loss <- loss_fn(outputs[[1]], batch_data, outputs[[2]], outputs[[3]])
          loss$backward()
          optimizer$step()
          total_loss <- total_loss + loss$item()
        }
      }
      
      denoised <- as.matrix(model(input)[[1]])
      
      list(model, denoised)
    },
  fun.extract = function(model)
    model[[2]],
  fun.apply_model = function(model, input) {
    as.matrix(model(input)[[1]])
  }
)