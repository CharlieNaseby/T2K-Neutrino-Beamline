import numpy as np
import matplotlib.pyplot as plt
import GPy as gp
from IPython.display import display




def drift(l):
    return np.array([[1, -2*l, l**2], [0, 1, -l], [0, 0, 1]], dtype=np.float64)


def focus_quad(sqrtk1, l):
    c = np.cos(l*sqrtk1)
    s = np.sin(l*sqrtk1)/sqrtk1

    cp = -sqrtk1*np.sin(l*sqrtk1)
    sp = np.cos(l*sqrtk1)
    return np.array([[c**2, -2*s*c, s**2], [-c*cp, s*cp + c*sp, -s*sp], [cp**2, -2*sp*cp, sp**2]], dtype = np.float64)


def defocus_quad(sqrtk1, l):
    c = np.cosh(l*sqrtk1)
    s = np.sinh(l*sqrtk1)/sqrtk1

    cp = sqrtk1*np.sinh(l*sqrtk1)
    sp = np.cosh(l*sqrtk1)
    return np.array([[c**2, -2*s*c, s**2], [-c*cp, s*cp + c*sp, -s*sp], [cp**2, -2*sp*cp, sp**2]], dtype = np.float64)
   
def get_gamma(alpha, beta):
    return (1+alpha**2)/beta

def get_width(emittance, beampars):
    return (emittance*beampars[0])**0.5

def get_acquisition_function(y, sigma, best_point, xi):
    return sigma*xi + (best_point-y)
    return sigma-(y-(best_point + xi))/(sigma+0.001)

def get_true_width(initial_beam, emittance, x):
    bl = Beamline()
    if(type(x) == np.float64):
        x = [x]
    resulting_beam = np.array([bl.propagate_beam(initial_beam, p) for p in x])
    return np.array([get_width(emittance, beam) for beam in resulting_beam])


class Beamline:
    strengths = []
    def __init__(self):
        self.strengths = 1.0
    
    def propagate_beam(self, initial_beam, strengths):
        output_beam = drift(1.0).dot(focus_quad(strengths, 1.0).dot(drift(2.0).dot(initial_beam)))
        return output_beam



def get_gp_func(variance, lengthscale, data_x, data_y):

    # Create RBF kernel with user-provided hyperparameters
    kernel = gp.kern.RBF(input_dim=1, variance=variance, lengthscale=lengthscale)
    
    # Create Gaussian Process model
    model = gp.models.GPRegression(data_x.reshape(-1, 1), data_y.reshape(-1, 1), kernel, noise_var=1e-10)
    
    model.Gaussian_noise.variance.fix(1e-10)
    model.optimize(messages=True)
    return model









emittance = 1
alpha = 0.1
beta = 1
gamma = get_gamma(alpha, beta)


initial_beam = np.array([beta, alpha, gamma])

predictor = Beamline()
final_beam = predictor.propagate_beam(initial_beam, 0.7)
target_width = get_width(emittance, final_beam)

print(f"target width = {target_width}")


x_train = np.linspace(0.01, 2, 4)
y_train = get_true_width(initial_beam, emittance, x_train)

x_truth = np.linspace(0.01, 2, 1000)
y_truth = get_true_width(initial_beam, emittance, x_truth)


    #hyperparmeters
variance = [16]
lengthscale = [0.7]
xi = 0.1

for i in range(20):
    squared_diff_train = (y_train-target_width)**2
    model = get_gp_func(variance, lengthscale, x_train, y_train)

    y_model_mean, y_model_var = model.predict(x_truth.reshape(-1, 1))
    quadratic_prediction = (y_model_mean-target_width)**2
    plus_onesig = (y_model_mean+np.sqrt(y_model_var)-target_width)**2
    best_point = np.argmin(quadratic_prediction)
    improvement = get_acquisition_function(quadratic_prediction, np.abs(plus_onesig - quadratic_prediction), np.min(quadratic_prediction), xi)


    plt.figure(figsize=(10, 6))
    plt.plot(x_truth, y_truth, 'b', label='Ground Truth')
    plt.plot(x_train, y_train, 'bo', label='Model Input')
    plt.plot(x_truth, y_model_mean, 'r-', label='Mean Prediction')
    plt.fill_between(x_truth.flatten(), (y_model_mean - np.sqrt(y_model_var)).flatten(), (y_model_mean + np.sqrt(y_model_var)).flatten(), color='pink', alpha=0.5, label='68% Credible Interval')
    plt.plot(x_truth, improvement, 'g', label='Acquisition Function')
    plt.plot(x_truth, (y_model_mean-target_width)**2, label='Squared Difference between prediction and target')
    plt.plot(x_truth, (y_model_mean.flatten()-y_truth), 'b', label='Model - True')
    plt.fill_between(x_truth.flatten(), (y_model_mean.flatten()-y_truth - np.sqrt(y_model_var).flatten()).flatten(), (y_model_mean.flatten()-y_truth + np.sqrt(y_model_var).flatten()).flatten(), color='blue', alpha=0.5)
    #plt.fill_between(x_truth.flatten(), (y_model_mean - np.sqrt(y_model_var)).flatten(), (y_model_mean + np.sqrt(y_model_var)).flatten(), color='pink', alpha=0.5, label='68% Credible Interval')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Gaussian Process Regression with RBF Kernel')
    plt.legend()
    plt.grid(True)
    plt.show()

    print(f"Best settings found at {x_truth[np.argmax(np.abs(y_model_mean-target_width))]}")

    x_train = np.append(x_train, x_truth[np.argmax(improvement)])
    y_train = np.append(y_train, get_true_width(initial_beam, emittance, x_truth[np.argmax(improvement)]))
    xi*=2 #increase the exploration as time goes on
    if(xi>10):
        xi=10




