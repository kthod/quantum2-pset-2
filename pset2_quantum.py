import matplotlib.pyplot as plt
import imageio # lib to make gifs
from scipy.special import hermite
from scipy.constants import pi
import numpy as np
from math import factorial
# %matplotlib inline
# #If your screen has retina display this will increase resolution of plots
# %config InlineBackend.figure_format = 'retina'

# Constants
hbar = 1
m = 1
omega = 1



# Time steps and time range
time_steps = 1000
t_max = 10
time_range = np.linspace(0, t_max, time_steps)

# Position range
x_range = np.linspace(-10, 10, 10000)

def N(v):
    '''Normalization constant '''
    
    return 1./np.sqrt(np.sqrt(np.pi)*2**v*factorial(v))


def psi(v, x):
    """Harmonic oscillator wavefunction for level v computed on grid of points x"""
    
    Hr=hermite(v)
    Psix = N(v)*Hr(x)*np.exp(-0.5*x**2)
    
    return Psix

# Store frames for the GIF
frames = []


# Initial state coefficients
c0 = 1 / np.sqrt(2)
c1 = 1 #/ np.sqrt(2)
c2 = 1 / np.sqrt(2)


# Calculate energy eigenvalues

# eigenvalues of |0> and |1>
E0 = 0.5 * hbar * omega
E1 = 1.5 * hbar * omega
E2 = 2.5 * hbar * omega
x_mean = []
x_mean_lectures = []
# Time evolution of the wavefunction
for t in time_range:
    
    exp_factor0 = np.exp(-1j * E0 * t / hbar)
    exp_factor1 = np.exp(-1j * E1 * t / hbar)
    exp_factor2 = np.exp(-1j * E2 * t / hbar)

    ct0 = c0 * exp_factor0
    ct1 = c1 * exp_factor1
    ct2 = c2 * exp_factor2

    psi_t =  ct1 * psi(1, x_range) #+ct2 * psi(2, x_range)
    prob_distribution = np.abs(psi_t)**2
    
    # Question : What is the meaning of the 2 following lines of code?
    
    # Calculate the mean value of the position operator
    x_t_mean = np.trapz(x_range * prob_distribution, x_range)
    x_mean.append(x_t_mean)
    
    # Calculate the mean value of the position operator - lectures
    x_t_mean_lectures = np.real(np.conj(psi_t) * x_range * psi_t)
    x_mean_lectures.append(x_t_mean)
    
    
    
    # # Plot probability distribution and save the frame
    # fig, ax = plt.subplots(figsize=(10, 5))
    # ax.plot(x_range, prob_distribution)
    # ax.set_xlabel('Position (x)')
    # ax.set_ylabel('Probability Distribution (|ψ(x,t)|²)')
    # ax.set_title(f'Time Evolution of the Probability Distribution for Quantum Harmonic Oscillator (t = {t:.2f})')

    # # Convert plot to image and store it
    # fig.canvas.draw()
    # image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    # image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    # frames.append(image)

    # # Close the plot to prevent display
    # plt.close(fig)

# Create a GIF using the frames
# GIF is saved to the current jupyter directory
# imageio.mimsave('quantum_harmonic_oscillator_superpos0_1.gif', frames, duration=5)



# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Plot x_mean on the first subplot
ax1.plot(time_range, x_mean)
ax1.set_xlabel('Time')
ax1.set_ylabel('Mean value of x')
ax1.set_title('Time Evolution of <x>')

# Plot x_mean_lectures on the second subplot
ax2.plot(time_range, x_mean_lectures)
ax2.set_xlabel('Time')
ax2.set_ylabel('Mean value of x')
ax2.set_title('Time Evolution of <x> (Lectures)')

# Adjust spacing between subplots
plt.subplots_adjust(wspace=0.3)

# Show the plot
plt.show()