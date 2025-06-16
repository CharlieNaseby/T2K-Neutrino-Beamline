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
from iminuit import Minuit

def get_gp_func(X, Y):

    likelihood = GaussianLikelihood()#noise_constraint=gpytorch.constraints.LessThan(1e-4))
    #likelihood.noise = 0.99e-4

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
NDIM_IN = 5  #how many free parameters
NDIM_OUT = 1 #the number of outputs the GP should predict, fixed to 1 for this approach


if __name__=='__main__':

    nseedpoints = 4
    bounds = np.array([[0.01, -2., 0.01, 0.01, 0.01],[2., -0.01, 1.0, 2, 1.0]])
    learn_rate = 1.0
    target_parameters = [0.7, -0.7, 0.1, 0.7, 0.2]
    target_properties = [[5.0, 3.0, 2.0, 2.0]] #target position x, position y, width x, width y 
    inter = Interface.Interface(target_parameters, target_properties)
    inter.SetInitialValues(target_parameters)
    inter.SetFileWriting(False)
    chisq = inter.MAE(target_parameters)
    target_properties = inter.GetBeamProperties()
    print(f'target properties {target_properties}')
#    inter.SetData(target_properties)

##### MINUIT section#####

#    minimiser_start_point = np.array([1., -1., 1.])
#    mini = Minuit(inter.fcn, minimiser_start_point)
#    mini.errordef = Minuit.LEAST_SQUARES
#    mini.limits = [(0.01, 2.0), (-2.0, -0.01), (0.01, 2.0)]
#
#    mini.migrad()
#    mini.hesse()
#
#    if(not mini.valid):
#        print(f'error minimum found is invalid')
#    print(f'minuit found a minimum at {mini.values} with error {mini.errors}')
#    minuit_fcn_history = inter.GetFcnHistory()


##### END MINUIT section#####

    #now set up for an asimov fit
    print(f'Set target beam properties to {target_parameters}')

    x_train = []
    for k in range(bounds[0].size):
        x_train.append(np.random.uniform(bounds[0,k], bounds[1,k], nseedpoints))

    x_train = np.array(x_train).T

#    x0_train = np.random.uniform(0.01, 2, nseedpoints)
#    x1_train = np.random.uniform(-2, -0.01, nseedpoints)
    
#    x_train = np.column_stack((x0_train, x1_train))
    y_train = np.array([-inter.MAE(x) for x in x_train])


    print(f'running simulations for strengths {x_train}')
    print(f'giving likelihood of {y_train}')
    
#    ntruepoints = 5
#
#    x_truth = np.linspace(0.01, 2, ntruepoints)
#    x0_truth, x1_truth = np.meshgrid(x_truth, -x_truth)
#    x_truth = np.column_stack((x0_truth.ravel(), x1_truth.ravel()))
#    y_truth = np.array([-inter.MAE(x) for x in x_truth])

#    solution_curve = []
#    for i in range(len(y_truth)):
#        if np.abs(y_truth[i])<1000:
#            solution_curve.append(x_truth[i])
#    
#    solution_curve = np.array(solution_curve)

    step = []
    best_chisq = []



    for i in range(100):
    
        with gpytorch.settings.cholesky_max_tries(10):
            try:
                model, mll = get_gp_func(torch.from_numpy(x_train), torch.from_numpy(-np.log10(-y_train)).reshape(-1, 1)) #just train x
            except:
                plt.plot(step, best_chisq, label='BO')
#                plt.plot(range(len(minuit_fcn_history)), np.log10(minuit_fcn_history), label='minuit')
                plt.xlabel('Iteration')
                plt.ylabel('MAE')
                plt.yscale('log')
                plt.legend()
                plt.show()
                print(f'final beam properties {inter.GetBeamProperties()}')
                break


  
    
    
        UCB = UpperConfidenceBound(model, beta=learn_rate)
        maximum = UpperConfidenceBound(model, beta=0.0)
    
    
        candidate, acq_value = optimize_acqf(
            UCB, bounds=torch.from_numpy(bounds), q=1, num_restarts=10, raw_samples=100,
        )

        best_fit, best_acq_value = optimize_acqf(
            maximum, bounds=torch.from_numpy(bounds), q=1, num_restarts=20, raw_samples=100,
        )

        print(f'Best fit parameters {best_fit} with acq_value {-best_acq_value}')
        with torch.no_grad():  # Ensure no gradients are computed
            posterior = model.posterior(best_fit.reshape(-1, NDIM_IN).to(torch.double))
    
        best_fit_chisq = posterior.mean
        best_fit_error = posterior.variance.sqrt()

        step.append(i)
        best_chisq.append(inter.MAE(best_fit.numpy().flatten()))

        print(f'With MAE pred = {-best_fit_chisq} +- {best_fit_error}')

        print(f'True MAE at this point = {best_chisq[-1]}')


#        x_array = torch.from_numpy(x_truth)
#        x_array_reshaped = x_array.reshape(-1, NDIM_IN).to(torch.double)
#    
#        with torch.no_grad():  # Ensure no gradients are computed
#            posterior = model.posterior(x_array_reshaped)
#    
#        mean = posterior.mean
#        ci_band = posterior.variance.sqrt()
#    
#    
#        x_array_reshaped = x_array.reshape(-1, 1, NDIM_IN).to(torch.double)
#        acquisition_values = UCB(x_array_reshaped)




#        print(f'Best proposal position = {x_truth[np.argmax(mean)]} with mean {mean[np.argmax(mean)]} and stddev {ci_band[np.argmax(mean)]}')
#        if(i%4==0):
#
#            fig, axes = plt.subplots(2, 3, figsize=(12, 5))
#
#            cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["red","white","blue"])
# 
#            norm=plt.Normalize(-4, np.log(-np.min(y_truth)))
#            true_surf = axes[0,0].pcolormesh(x0_truth, x1_truth, np.log(-y_truth.reshape(ntruepoints,ntruepoints)), cmap='viridis', norm=norm)
#            plt.colorbar(true_surf, label='Log(True Chisq)')
#            axes[0,0].scatter(target_parameters[0], target_parameters[1], c='black', label='True Target', s=0.3)
#            #axes[0,0].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
#            axes[0,0].set_title('True Chisq')
#            axes[0,0].set_ylabel('Q2 k1')
#        
#            gp_surf = axes[0,1].pcolormesh(x0_truth, x1_truth, np.log(np.abs(mean.reshape(ntruepoints,ntruepoints))), cmap='viridis', norm=norm)
#            axes[0,1].scatter(x_train[:,0], x_train[:,1], color='red', s=0.3)
#        #    axes[0,1].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
#            fig.colorbar(gp_surf, label='Log(abs(Model Beam Chisq))')
#            axes[0,1].set_title('Chisq Predicted from GP')
#        
#        
#
##            norm=plt.Normalize(0.0 ,np.log(-np.max(np.abs(y_truth-mean.numpy().flatten()))))
#            diff_surf = axes[0,2].pcolormesh(x0_truth, x1_truth, np.log(np.abs((y_truth-mean.numpy().flatten()).reshape(ntruepoints,ntruepoints))), cmap='viridis', norm=norm)
#            axes[0,2].scatter(x_train[:,0], x_train[:,1], color='red', s=0.3)
#        #    axes[0,2].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
#            fig.colorbar(diff_surf, label='Log(abs(Chisq True-Model))')
#            axes[0,2].set_title('True - Prediction Chisq')
#
#
#            norm=plt.Normalize(0.0,np.max(np.abs(ci_band.numpy().flatten())))
#            uncert_surf = axes[1,0].pcolormesh(x0_truth, x1_truth, (ci_band.numpy().flatten()).reshape(ntruepoints,ntruepoints), cmap='viridis', norm=norm)
#            axes[1,0].scatter(x_train[:,0], x_train[:,1], color='red', s=0.3)
#        #    axes[1,0].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
#            fig.colorbar(uncert_surf, label='Beam Model Uncert')
#            axes[1,0].set_title('Uncertainty')
#            axes[1,0].set_xlabel('Q1 k1')
#            axes[1,0].set_ylabel('Q2 k1')
#
#
#
#            norm=plt.Normalize(np.min(acquisition_values.detach().numpy().flatten()) ,np.max(acquisition_values.detach().numpy().flatten()))
#            acquisition_surf = axes[1,1].pcolormesh(x0_truth, x1_truth, ((np.abs(acquisition_values.detach().numpy())).flatten()).reshape(ntruepoints,ntruepoints), cmap='viridis', norm=norm)
#            axes[1,1].scatter(x_train[:,0], x_train[:,1], color='red', s=0.3)
#            axes[1,1].scatter(x_truth[np.argmax(acquisition_values.detach().numpy()),0], x_truth[np.argmax(acquisition_values.detach().numpy()),1])
#        #    axes[1,1].scatter(solution_curve[:,0], solution_curve[:,1], c='blue', s=0.1)
#            fig.colorbar(acquisition_surf, label='Acquisition Function')
#            axes[1,1].set_title('Acquisition Function')
#
#            axes[1,2].plot(step, np.log(best_chisq))
#            axes[1,2].set_xlabel('Iteration')
#            axes[1,2].set_ylabel('Log(Chisq)')
#
#            plt.tight_layout()
#            plt.show()

        x_train = np.append(x_train, candidate.numpy(), axis=0)
        y_train = np.append(y_train, np.array([-inter.MAE(x_train[-1])]), axis=0)
        if(y_train.size > 40):
            x_train = x_train[-40:]
            y_train = y_train[-40:]
        print(f'y train is now {y_train}')
        learn_rate = learn_rate*0.99

#    if(not mini.valid):
#        print(f'error minimum found is invalid')
#    print(f'minuit found a minimum at {mini.values} with error {mini.errors} and fcn {inter.fcn(mini.values)}')

    plt.plot(step, best_chisq, label='BO')
#    plt.plot(range(len(minuit_fcn_history)), minuit_fcn_history, label='minuit')
    plt.xlabel('Iteration')
    plt.ylabel('MAE')
    plt.yscale('log')
    plt.legend()
    plt.show()


