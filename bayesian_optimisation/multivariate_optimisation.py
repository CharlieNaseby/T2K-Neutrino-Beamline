import torch
from torch import Tensor
import gpytorch
from gpytorch.likelihoods import GaussianLikelihood
from botorch.models import SingleTaskGP
from botorch.models.transforms import Normalize, Standardize
from botorch.fit import fit_gpytorch_mll
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.acquisition import LogExpectedImprovement
from botorch.optim import optimize_acqf
from botorch.acquisition import AcquisitionFunction
from botorch.models.model import Model
from botorch.utils.transforms import t_batch_mode_transform
from botorch.acquisition.monte_carlo import MCAcquisitionFunction
from botorch.sampling.base import MCSampler
from botorch.sampling.normal import SobolQMCNormalSampler

import matplotlib.pyplot as plt
import numpy as np
import GPy as gp
from IPython.display import display
import math
from typing import Optional

def get_transfer_mtx(c, s, cp, sp):
    return np.array([[c**2, -2*s*c, s**2], [-c*cp, s*cp + c*sp, -s*sp], [cp**2, -2*sp*cp, sp**2]], dtype = np.float64)

def drift(l):
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
   
def get_gamma(alpha, beta):
    return (1+alpha**2)/beta

def get_width(emittance, beampars):
    return (emittance*beampars[0,:])**0.5

def get_acquisition_function(y, sigma, best_point, xi):
    return sigma*xi + (best_point-y)

def get_true_width(initial_beam, emittance, x):
    bl = Beamline()
    if(type(x) == np.float64 or x.shape == (NDIM_IN,)):
        x = [x]
    resulting_beam = np.array([bl.propagate_beam(initial_beam, p) for p in x])
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
      input_transform=Normalize(d=NDIM_IN),
      outcome_transform=Standardize(m=NDIM_OUT),
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

        square_nsigma_to_target = ((mean-target_width[0])/sigma)**2 ##TODO CERN
        sigma_penalty = learning_rate/sigma**2
        return -(square_nsigma_to_target + sigma_penalty).squeeze()


torch.set_default_dtype(torch.float64)
np.random.seed(1989)


NDIM_IN = 2
NDIM_OUT = 1

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



for i in range(10):
    model, mll = get_gp_func(torch.from_numpy(x_train), torch.from_numpy(y_train[:,0]).reshape(-1, 1)) #just train x

    custom_log_ei = qScalarizedUpperConfidenceBound(model=model, beta=0.1, weights=torch.tensor([0.1]))


    bounds = torch.tensor([[0.01, -2.],[2., -0.01]])

    custom_candidate, custom_acq_value = optimize_acqf(
      custom_log_ei, bounds=bounds, q=1, num_restarts=5, raw_samples=20,
    )

    x_array = torch.from_numpy(x_truth)
    x_array_reshaped = x_array.reshape(-1, NDIM_IN).to(torch.double)

    with torch.no_grad():  # Ensure no gradients are computed
        posterior = model.posterior(x_array_reshaped)

    mean = posterior.mean
    ci_band = posterior.variance.sqrt()


    x_array_reshaped = x_array.reshape(-1, 1, NDIM_IN).to(torch.double)
    custom_acquisition_values = custom_log_ei(x_array_reshaped)

    fig, axes = plt.subplots(2, 3, figsize=(12, 5))

    widthmap_x = axes[0,0].scatter(x_truth[:,0], x_truth[:,1], c=y_truth[:,0], cmap='viridis',s=0.3)
    fig.colorbar(widthmap_x, label='True Beam Width Arb.')
    axes[0,0].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
    axes[0,0].set_title('X width')
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


    widthmap_x_gp = axes[0,2].scatter(x_truth[:,0], x_truth[:,1], c=y_truth[:,0]-mean.numpy()[:,0], cmap='viridis',s=0.3)
    axes[0,2].scatter(x_train[:,0], x_train[:,1], color='red', s=0.1)
    axes[0,2].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
    fig.colorbar(widthmap_x_gp, label='Beam Width True-Model')
    axes[0,2].set_title('True - Prediction X')

    widthmap_xsigma_gp = axes[1,0].scatter(x_truth[:,0], x_truth[:,1], c=ci_band.numpy()[:,0], cmap='viridis', s=0.3)
    axes[1,0].scatter(x_train[:,0], x_train[:,1], color='red', s=0.1)
    axes[1,0].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
    fig.colorbar(widthmap_xsigma_gp, label='Beam Model Uncert')
    axes[1,0].set_title('Uncertainty')
    axes[1,0].set_xlabel('Q1 k1**0.5')
    axes[1,0].set_ylabel('Q2 k1**0.5')

    acquisition = axes[1,1].scatter(x_truth[:,0], x_truth[:,1], c=-np.log(np.abs(custom_acquisition_values.detach().numpy())), cmap='viridis', s=0.3)
    axes[1,1].scatter(x_train[:,0], x_train[:,1], color='red', s=0.1)
    axes[1,1].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
    fig.colorbar(acquisition, label='log(Acquisition Function)')
    axes[1,1].set_title('Acquisition Function')
    plt.tight_layout()
    plt.show()


    x_train = np.append(x_train, custom_candidate.numpy(), axis=0)
    y_train = np.append(y_train, [get_true_width(initial_beam, emittance, x_train[-1])[0]], axis=0)
