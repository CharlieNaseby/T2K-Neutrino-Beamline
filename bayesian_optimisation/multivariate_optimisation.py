import torch
import gpytorch
from botorch.models import SingleTaskGP
from botorch.models.transforms import Normalize, Standardize
from botorch.fit import fit_gpytorch_mll
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.acquisition import LogExpectedImprovement
from botorch.optim import optimize_acqf
import matplotlib.pyplot as plt
import numpy as np
from gpytorch.likelihoods import GaussianLikelihood
from botorch.acquisition import AcquisitionFunction
from botorch.models.model import Model
from botorch.utils.transforms import t_batch_mode_transform
from botorch.acquisition.monte_carlo import MCAcquisitionFunction
from torch import Tensor
from botorch.sampling.base import MCSampler
from botorch.sampling.normal import SobolQMCNormalSampler
import GPy as gp
from IPython.display import display
import math
from typing import Optional

def get_transfer_mtx(c, s, cp, sp):
    return np.array([[c**2, -2*s*c, s**2], [-c*cp, s*cp + c*sp, -s*sp], [cp**2, -2*sp*cp, sp**2]], dtype = np.float64)

def drift(l):
#    return np.array([[1, -2*l, l**2], [0, 1, -l], [0, 0, 1]], dtype=np.float64)
    return get_transfer_mtx(1, l, 0, 1)


def focus_quad(sqrtk1, l):
    c = np.cos(l*sqrtk1)
    s = np.sin(l*sqrtk1)/sqrtk1

    cp = -sqrtk1*np.sin(l*sqrtk1)
    sp = np.cos(l*sqrtk1)

    ch = np.cos(l*sqrtk1)
    sh = np.sin(l*sqrtk1)/sqrtk1

    chp = sqrtk1*np.sin(l*sqrtk1)
    shp = np.cos(l*sqrtk1)
    return np.array([get_transfer_mtx(c, s, cp, sp), get_transfer_mtx(ch, sh, chp, shp)])
#    return np.array([[c**2, -2*s*c, s**2], [-c*cp, s*cp + c*sp, -s*sp], [cp**2, -2*sp*cp, sp**2]], dtype = np.float64)


def defocus_quad(sqrtk1, l):
    c = np.cos(l*sqrtk1)
    s = np.sin(l*sqrtk1)/sqrtk1

    cp = -sqrtk1*np.sin(l*sqrtk1)
    sp = np.cos(l*sqrtk1)

    ch = np.cos(l*sqrtk1)
    sh = np.sin(l*sqrtk1)/sqrtk1

    chp = sqrtk1*np.sin(l*sqrtk1)
    shp = np.cos(l*sqrtk1)

    return np.array([get_transfer_mtx(ch, sh, chp, shp), get_transfer_mtx(c, s, cp, sp)])
#    return np.array([[c**2, -2*s*c, s**2], [-c*cp, s*cp + c*sp, -s*sp], [cp**2, -2*sp*cp, sp**2]], dtype = np.float64)
   
def get_gamma(alpha, beta):
    return (1+alpha**2)/beta

def get_width(emittance, beampars):
    return (emittance*beampars[0,:])**0.5

def get_acquisition_function(y, sigma, best_point, xi):
    return sigma*xi + (best_point-y)
    return sigma-(y-(best_point + xi))/(sigma+0.001)

def get_true_width(initial_beam, emittance, x):
    bl = Beamline()
    if(type(x) == np.float64 or x.shape == (2,)):
        x = [x]
    resulting_beam = np.array([bl.propagate_beam(initial_beam, p) for p in x])
#    print(f'resulting beam = {resulting_beam}')
    return np.array([get_width(emittance, beam) for beam in resulting_beam])


class Beamline:
    strengths = []
    def __init__(self):
        self.strengths = 1.0
    def apply_quad(self, strength, l, beam):
        if(strength<0.0):
            quad = defocus_quad(-strength, l)
        else:
            quad = focus_quad(strength, l)
        return np.array([quad[0].dot(beam[:,0]), quad[1].dot(beam[:,1])]).T

    def propagate_beam(self, initial_beam, strengths):
        output_beam = drift(1.0).dot(self.apply_quad(strengths[1], 1.0, drift(1.0).dot(self.apply_quad(strengths[0], 1.0, drift(1.0).dot(initial_beam)))))
        return output_beam



def get_gp_func(X, Y):

    likelihood = GaussianLikelihood(noise_constraint=gpytorch.constraints.LessThan(1e-4))
    likelihood.noise = 0.99e-4

    gp = SingleTaskGP(
      train_X=X,
      train_Y=Y,
      input_transform=Normalize(d=2),
      outcome_transform=Standardize(m=1),
      likelihood=likelihood)
    mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
    fit_gpytorch_mll(mll)

    return gp, mll

class qScalarizedUpperConfidenceBound(MCAcquisitionFunction):
    def __init__(
        self,
        model: Model,
        beta: Tensor,
        weights: Tensor,
        sampler: Optional[MCSampler] = None,
    ) -> None:
        # we use the AcquisitionFunction constructor, since that of
        # MCAcquisitionFunction performs some validity checks that we don't want here
        super(MCAcquisitionFunction, self).__init__(model=model)
        if sampler is None:
            sampler = SobolQMCNormalSampler(sample_shape=torch.Size([512]))
        self.sampler = sampler
        self.register_buffer("beta", torch.as_tensor(beta))
        self.register_buffer("weights", torch.as_tensor(weights))

    @t_batch_mode_transform()
    def forward(self, X: Tensor) -> Tensor:
        """Evaluate scalarized qUCB on the candidate set `X`.

        Args:
            X: A `(b) x q x d`-dim Tensor of `(b)` t-batches with `q` `d`-dim
                design points each.

        Returns:
            Tensor: A `(b)`-dim Tensor of Upper Confidence Bound values at the
                given design points `X`.
        """
        posterior = self.model.posterior(X)

        mean = posterior.mean
        sigma = posterior.variance.sqrt()
        learning_rate=0.1
#        target = torch.from_numpy(target_width[0])

#        print(f'mean shape = {mean.shape} target width shape = {target.shape}')
        square_nsigma_to_target = ((mean-target_width[0])/sigma)**2 ##TODO CERN
        sigma_penalty = learning_rate/sigma**2
        return -(square_nsigma_to_target + sigma_penalty).squeeze()
        print(f'shapes x = {X.shape}, sigma = {sigma.shape}')
        print(f'test return shape = {sigma.squeeze().shape}')
        print(sigma.squeeze())
        return sigma.squeeze()

        samples = self.get_posterior_samples(posterior)  # n x b x q x o
        scalarized_samples = samples.matmul(self.weights)  # n x b x q
        mean = posterior.mean  # b x q x o
        scalarized_mean = mean.matmul(self.weights)  # b x q
        ucb_samples = (
            scalarized_mean
            + math.sqrt(self.beta * math.pi / 2)
            * (scalarized_samples - scalarized_mean).abs()
        )
        print(f'working ret type {ucb_samples.max(dim=-1)[0].mean(dim=0).shape}')
        return ucb_samples.max(dim=-1)[0].mean(dim=0)




torch.set_default_dtype(torch.float64)


np.random.seed(1989)

def true_function(x1, x2):
    return np.sin(x1)# + np.cos(x2)

# Generate random input data
np.random.seed(0)
n = 50  # Number of data points
X1 = np.random.uniform(0, 2 * np.pi, n)
X2 = np.random.uniform(0, 2 * np.pi, n)
Y = true_function(X1, X2) + np.random.normal(0, 0.1, n)  # Add some noise

# Convert to PyTorch tensors
train_X = torch.tensor(np.column_stack((X1, X2)), dtype=torch.float64)
train_Y = torch.tensor(Y, dtype=torch.float64).reshape(-1, 1)
X1 = torch.tensor(train_X)#.reshape(-1, 2)

# Step 2: Fit a Gaussian Process model
model = SingleTaskGP(X1, train_Y)

# Use ExactMarginalLogLikelihood for fitting the model
mll = ExactMarginalLogLikelihood(model.likelihood, model)


bounds = torch.tensor([[0., 0],[2*np.pi, 2*np.pi]])
#bounds = torch.tensor([[0.],[ 2*np.pi]])

logEI = LogExpectedImprovement(model=model, best_f=train_Y.max())
candidate, acq_value = optimize_acqf(
  logEI, bounds=bounds, q=1, num_restarts=5, raw_samples=20,
)

# Fit the model
def fit_gpytorch_model(mll):
    # Use the Adam optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=0.1)

    # Training loop
    for _ in range(100):
        optimizer.zero_grad()
        output = model(train_X)
        loss = -mll(output, train_Y)
        loss.backward()
        optimizer.step()

#fit_gpytorch_model(mll)


emittance = np.array([1., 1.])
alpha = np.array([0.1, 0.1])
beta = np.array([1., 1.])
gamma = get_gamma(alpha, beta)


initial_beam = np.array([beta, alpha, gamma])


predictor = Beamline()
final_beam = predictor.propagate_beam(initial_beam, np.array([0.7, -0.7]))
target_width = get_width(emittance, final_beam)


x_train = np.linspace(0.01, 2, 4, dtype=np.float64)
x0_train, x1_train = np.meshgrid(x_train, -x_train)
x_train = np.column_stack((x0_train.ravel(), x1_train.ravel()))
y_train = get_true_width(initial_beam, emittance, x_train)



x_truth = np.linspace(0.01, 2, 100)
x0_truth, x1_truth = np.meshgrid(x_truth, -x_truth)
x_truth = np.column_stack((x0_truth.ravel(), x1_truth.ravel()))
y_truth = get_true_width(initial_beam, emittance, x_truth)

solution_curve = []
for i in range(len(y_truth)):
    if np.abs((y_truth[i]-target_width)[0])<0.05:
        solution_curve.append(x_truth[i])

solution_curve = np.array(solution_curve)

#fig, axes = plt.subplots(1, 2, figsize=(12, 5))
#
#widthmap_x = axes[0].scatter(x_truth[:,0], x_truth[:,1], c=y_truth[:,0], cmap='viridis')
#fig.colorbar(widthmap_x, label='Beam Width Arb.')
#axes[0].set_title('X width')
#axes[0].set_xlabel('Q1 k1**0.5')
#axes[0].set_ylabel('Q2 k1**0.5')
#
#widthmap_y = axes[1].scatter(x_truth[:,0], x_truth[:,1], c=y_truth[:,1], cmap='viridis')
#fig.colorbar(widthmap_y, label='Beam Width Arb.')
#axes[1].set_title('Y width')
#axes[1].set_xlabel('Q1 k1**0.5')
#axes[1].set_ylabel('Q2 k1**0.5')

#plt.show()


#plt.plot(x_truth, y_truth[:,0], label='x width')
#plt.plot(x_truth, y_truth[:,1], label='y width')
#plt.legend()
#plt.show()

for i in range(10):
    model, mll = get_gp_func(torch.from_numpy(x_train), torch.from_numpy(y_train[:,0]).reshape(-1, 1)) #just train x

    logEI = LogExpectedImprovement(model=model, best_f=y_train.max())
    custom_log_ei = qScalarizedUpperConfidenceBound(model=model, beta=0.1, weights=torch.tensor([0.1]))


    bounds = torch.tensor([[0.01, -2.],[2., -0.01]])
    #bounds = torch.stack([torch.zeros(1), 2.*torch.ones(1)]).to(torch.double)
#    candidate, acq_value = optimize_acqf(
#      logEI, bounds=bounds, q=1, num_restarts=5, raw_samples=20,
#    )

    custom_candidate, custom_acq_value = optimize_acqf(
      custom_log_ei, bounds=bounds, q=1, num_restarts=5, raw_samples=20,
    )


    x_array = torch.from_numpy(x_truth)
    x_array_reshaped = x_array.reshape(-1, 2).to(torch.double)
    # Compute the posterior for the input array
    with torch.no_grad():  # Ensure no gradients are computed
        posterior = model.posterior(x_array_reshaped)

    mean = posterior.mean
    ci_band = posterior.variance.sqrt()


    x_array_reshaped = x_array.reshape(-1, 1, 2).to(torch.double)
    acquisition_values = logEI(x_array_reshaped)
    custom_acquisition_values = custom_log_ei(x_array_reshaped)

    fig, axes = plt.subplots(2, 3, figsize=(12, 5))

    widthmap_x = axes[0,0].scatter(x_truth[:,0], x_truth[:,1], c=y_truth[:,0], cmap='viridis',s=0.3)
    fig.colorbar(widthmap_x, label='True Beam Width Arb.')
    axes[0,0].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
    axes[0,0].set_title('X width')
    axes[0,0].set_xlabel('Q1 k1**0.5')
    axes[0,0].set_ylabel('Q2 k1**0.5')

#    widthmap_y = axes[0,1].scatter(x_truth[:,0], x_truth[:,1], c=y_truth[:,1], cmap='viridis')
#    fig.colorbar(widthmap_y, label='Beam Width Arb.')
#    axes[0,1].set_title('Y width')
#    axes[0,1].set_xlabel('Q1 k1**0.5')
#    axes[0,1].set_ylabel('Q2 k1**0.5')

    widthmap_x_gp = axes[0,1].scatter(x_truth[:,0], x_truth[:,1], c=mean, cmap='viridis',s=0.3)
    axes[0,1].scatter(x_train[:,0], x_train[:,1], color='red', s=0.1)
    axes[0,1].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
    fig.colorbar(widthmap_x_gp, label='Model Beam Width Arb.')
    axes[0,1].set_title('X width Predicted from GP')
    axes[0,1].set_xlabel('Q1 k1**0.5')
    axes[0,1].set_ylabel('Q2 k1**0.5')


    widthmap_x_gp = axes[0,2].scatter(x_truth[:,0], x_truth[:,1], c=y_truth[:,0]-mean.numpy()[:,0], cmap='viridis',s=0.3)
    axes[0,2].scatter(x_train[:,0], x_train[:,1], color='red', s=0.1)
    axes[0,2].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
    fig.colorbar(widthmap_x_gp, label='Beam Width True-Model')
    axes[0,2].set_title('True - Prediction X')
    axes[0,2].set_xlabel('Q1 k1**0.5')
    axes[0,2].set_ylabel('Q2 k1**0.5')

    widthmap_xsigma_gp = axes[1,0].scatter(x_truth[:,0], x_truth[:,1], c=ci_band.numpy()[:,0], cmap='viridis', s=0.3)
    axes[1,0].scatter(x_train[:,0], x_train[:,1], color='red', s=0.1)
    axes[1,0].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
    fig.colorbar(widthmap_xsigma_gp, label='Beam Model Uncert')
    axes[1,0].set_title('Uncertainty')
    axes[1,0].set_xlabel('Q1 k1**0.5')
    axes[1,0].set_ylabel('Q2 k1**0.5')

    print(custom_acquisition_values.detach().numpy().shape)
    acquisition = axes[1,1].scatter(x_truth[:,0], x_truth[:,1], c=custom_acquisition_values.detach().numpy(), cmap='viridis', s=0.3)
    axes[1,1].scatter(x_train[:,0], x_train[:,1], color='red', s=0.1)
    axes[1,1].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
    fig.colorbar(acquisition, label='Acquisition Function')
    axes[1,1].set_title('Acquisition Function')
    axes[1,1].set_xlabel('Q1 k1**0.5')
    axes[1,1].set_ylabel('Q2 k1**0.5')
    plt.show()



#    plt.plot(x_truth, mean.flatten(), color='red', label='Mean Prediction')
#    plt.fill_between(x_truth.flatten(), (mean - ci_band).flatten(), (mean+ci_band).flatten(), color='pink', alpha=0.5, label='68% Credible Interval')
#    plt.scatter(x_train, y_train, label='Training Data')
#    plt.plot(x_truth, y_truth, label='Ground Truth')
##    plt.scatter(candidate[0], acq_value, label='Selected sample point')
#    plt.plot(x_truth, acquisition_values.detach().numpy(), label='Acquisition Function')
#    plt.plot(x_truth, target_width*np.ones(x_truth.shape), label='Target')
#    plt.plot(x_truth, custom_acquisition_values.detach().numpy(), label='Custom Acquisition Function')
#    plt.scatter(custom_candidate[0], custom_acq_value, label='Selected custom sample point')
#    plt.legend()
#    plt.show()

    x_train = np.append(x_train, custom_candidate.numpy(), axis=0)
    y_train = np.append(y_train, [get_true_width(initial_beam, emittance, x_train[-1])[0]], axis=0)


    #hyperparmeters
#variance = [16]
#lengthscale = [0.7]
#xi = 0.1
#
#for i in range(20):
#    squared_diff_train = (y_train-target_width)**2
#    model = get_gp_func(variance, lengthscale, x_train, y_train)
#
#    y_model_mean, y_model_var = model.predict(x_truth.reshape(-1, 1))
#    quadratic_prediction = (y_model_mean-target_width)**2
#    plus_onesig = (y_model_mean+np.sqrt(y_model_var)-target_width)**2
#    best_point = np.argmin(quadratic_prediction)
#    improvement = get_acquisition_function(quadratic_prediction, np.abs(plus_onesig - quadratic_prediction), np.min(quadratic_prediction), xi)
#
#
#    plt.figure(figsize=(10, 6))
#    plt.plot(x_truth, y_truth, 'b', label='Ground Truth')
#    plt.plot(x_train, y_train, 'bo', label='Model Input')
#    plt.plot(x_truth, y_model_mean, 'r-', label='Mean Prediction')
#    plt.fill_between(x_truth.flatten(), (y_model_mean - np.sqrt(y_model_var)).flatten(), (y_model_mean + np.sqrt(y_model_var)).flatten(), color='pink', alpha=0.5, label='68% Credible Interval')
#    plt.plot(x_truth, improvement, 'g', label='Acquisition Function')
#    plt.plot(x_truth, (y_model_mean-target_width)**2, label='Squared Difference between prediction and target')
#    plt.plot(x_truth, (y_model_mean.flatten()-y_truth), 'b', label='Model - True')
#    plt.fill_between(x_truth.flatten(), (y_model_mean.flatten()-y_truth - np.sqrt(y_model_var).flatten()).flatten(), (y_model_mean.flatten()-y_truth + np.sqrt(y_model_var).flatten()).flatten(), color='blue', alpha=0.5)
#    #plt.fill_between(x_truth.flatten(), (y_model_mean - np.sqrt(y_model_var)).flatten(), (y_model_mean + np.sqrt(y_model_var)).flatten(), color='pink', alpha=0.5, label='68% Credible Interval')
#    plt.xlabel('x')
#    plt.ylabel('y')
#    plt.title('Gaussian Process Regression with RBF Kernel')
#    plt.legend()
#    plt.grid(True)
#    plt.show()
#
#    print(f"Best settings found at {x_truth[np.argmax(np.abs(y_model_mean-target_width))]}")
#
#    x_train = np.append(x_train, x_truth[np.argmax(improvement)])
#    y_train = np.append(y_train, get_true_width(initial_beam, emittance, x_truth[np.argmax(improvement)]))
#    xi*=2 #increase the exploration as time goes on
#    if(xi>10):
#        xi=10




