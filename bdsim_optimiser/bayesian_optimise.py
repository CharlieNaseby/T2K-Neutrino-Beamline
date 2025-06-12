import Interface
import numpy as np
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
from botorch.acquisition import UpperConfidenceBound
from botorch.settings import debug

import matplotlib.pyplot as plt
import matplotlib as mpl
import GPy as gp
from IPython.display import display
import math
from typing import Optional
from multiprocessing import Pool

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
np.random.seed(1987)

#global constants
NDIM_IN = 2  #how many free parameters
NDIM_OUT = 1 #the number of outputs the GP should predict, fixed to 1 for this approach
NSSEM = 2 #number of ssem planes to use 
NAXIS = 2 #number of spatial axes to sum over in the likelihood max of 2 for x and y


if __name__=='__main__':


    target_parameters = [0.7, -0.7]
    target_properties = [[0, 0, 1.0, 1.0]] #target position x, position y, width x, width y
    inter = Interface.Interface(target_parameters, target_properties)
    inter.SetInitialValues(target_parameters)
    inter.SetFileWriting(False)
    chisq = inter.fcn(target_parameters)
    target_properties = inter.GetBeamProperties()
    inter.SetData(target_properties)

    #now set up for an asimov fit
    print(f'Set target beam properties to {target_parameters}')

    x0_train = np.random.uniform(0.01, 2, 4)
    x1_train = np.random.uniform(-2, -0.01, 4)
    
    x_train = np.column_stack((x0_train, x1_train))
    y_train = np.array([-inter.fcn(x) for x in x_train])


    print(f'running simulations for strengths {x_train}')
    print(f'giving likelihood of {y_train}')
    
    ntruepoints = 20

    x_truth = np.linspace(0.01, 2, ntruepoints)
    x0_truth, x1_truth = np.meshgrid(x_truth, -x_truth)
    x_truth = np.column_stack((x0_truth.ravel(), x1_truth.ravel()))
    y_truth = np.array([-inter.fcn(x) for x in x_truth])

    solution_curve = []
    for i in range(len(y_truth)):
        if np.abs(y_truth[i])<1000:
            solution_curve.append(x_truth[i])
    
    solution_curve = np.array(solution_curve)
    for i in range(40):
    
        with gpytorch.settings.cholesky_max_tries(10):
            model, mll = get_gp_func(torch.from_numpy(x_train), torch.from_numpy(y_train).reshape(-1, 1)) #just train x
    
    
        UCB = UpperConfidenceBound(model, beta=2.0)
        logEI = LogExpectedImprovement(model=model, best_f=y_train.max())
        custom_log_ei = qScalarizedUpperConfidenceBound(model=model, beta=0.1, weights=torch.tensor([0.1]))
    
    
        bounds = torch.tensor([[0.01, -2.],[2., -0.01]])
    
        candidate, acq_value = optimize_acqf(
          logEI, bounds=bounds, q=1, num_restarts=25, raw_samples=200,
        )
    
        candidate, acq_value = optimize_acqf(
            UCB, bounds=bounds, q=1, num_restarts=10, raw_samples=100,
        )
    
        x_array = torch.from_numpy(x_truth)
        x_array_reshaped = x_array.reshape(-1, NDIM_IN).to(torch.double)
    
        with torch.no_grad():  # Ensure no gradients are computed
            posterior = model.posterior(x_array_reshaped)
    
        mean = posterior.mean
        ci_band = posterior.variance.sqrt()
    
    
        x_array_reshaped = x_array.reshape(-1, 1, NDIM_IN).to(torch.double)
    
        acquisition_values = logEI(x_array_reshaped)
    
        acquisition_values = UCB(x_array_reshaped)


        print(f'Best found position = {x_truth[np.argmax(mean)]} with mean {mean[np.argmax(mean)]} and stddev {ci_band[np.argmax(mean)]}')
        if(i%4==0):

            fig, axes = plt.subplots(2, 3, figsize=(12, 5))

            cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["red","white","blue"])
 
            norm=plt.Normalize(-4, np.log(-np.min(y_truth)))
            true_surf = axes[0,0].pcolormesh(x0_truth, x1_truth, np.log(-y_truth.reshape(ntruepoints,ntruepoints)), cmap='viridis', norm=norm)
            plt.colorbar(true_surf, label='Log(True Chisq)')
            axes[0,0].scatter(target_parameters[0], target_parameters[1], c='black', label='True Target', s=0.3)
            #axes[0,0].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
            axes[0,0].set_title('True Chisq')
            axes[0,0].set_ylabel('Q2 k1')
        
            gp_surf = axes[0,1].pcolormesh(x0_truth, x1_truth, np.log(np.abs(mean.reshape(ntruepoints,ntruepoints))), cmap='viridis', norm=norm)
            axes[0,1].scatter(x_train[:,0], x_train[:,1], color='red', s=0.3)
        #    axes[0,1].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
            fig.colorbar(gp_surf, label='Log(abs(Model Beam Chisq))')
            axes[0,1].set_title('Chisq Predicted from GP')
        
        

#            norm=plt.Normalize(0.0 ,np.log(-np.max(np.abs(y_truth-mean.numpy().flatten()))))
            diff_surf = axes[0,2].pcolormesh(x0_truth, x1_truth, np.log(np.abs((y_truth-mean.numpy().flatten()).reshape(ntruepoints,ntruepoints))), cmap='viridis', norm=norm)
            axes[0,2].scatter(x_train[:,0], x_train[:,1], color='red', s=0.3)
        #    axes[0,2].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
            fig.colorbar(diff_surf, label='Log(abs(Chisq True-Model))')
            axes[0,2].set_title('True - Prediction Chisq')


            norm=plt.Normalize(0.0,np.max(np.abs(ci_band.numpy().flatten())))
            uncert_surf = axes[1,0].pcolormesh(x0_truth, x1_truth, (ci_band.numpy().flatten()).reshape(ntruepoints,ntruepoints), cmap='viridis', norm=norm)
            axes[1,0].scatter(x_train[:,0], x_train[:,1], color='red', s=0.3)
        #    axes[1,0].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
            fig.colorbar(uncert_surf, label='Beam Model Uncert')
            axes[1,0].set_title('Uncertainty')
            axes[1,0].set_xlabel('Q1 k1')
            axes[1,0].set_ylabel('Q2 k1')



            norm=plt.Normalize(np.min(acquisition_values.detach().numpy().flatten()) ,np.max(acquisition_values.detach().numpy().flatten()))
            acquisition_surf = axes[1,1].pcolormesh(x0_truth, x1_truth, ((np.abs(acquisition_values.detach().numpy())).flatten()).reshape(ntruepoints,ntruepoints), cmap='viridis', norm=norm)
            axes[1,1].scatter(x_train[:,0], x_train[:,1], color='red', s=0.3)
            axes[1,1].scatter(x_truth[np.argmax(acquisition_values.detach().numpy()),0], x_truth[np.argmax(acquisition_values.detach().numpy()),1])
        #    axes[1,1].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
            fig.colorbar(acquisition_surf, label='Acquisition Function')
            axes[1,1].set_title('Acquisition Function')
            plt.tight_layout()
            plt.show()

        x_train = np.append(x_train, candidate.numpy(), axis=0)
        y_train = np.append(y_train, np.array([-inter.fcn(x_train[-1])]), axis=0)
        print(f'y train is now {y_train}')

