import numpy as np
import matplotlib.pyplot as plt

# Given data: node counts and corresponding L2 errors
node_counts = np.array([6, 15, 45])
errors = np.array([1.278357e-02, 6.140789e-03, 2.240537e-03])

# Estimate mesh size h ~ 1/sqrt(N)
h_values = 1 / np.sqrt(node_counts)

# Take log10 for both h and error
log_h = np.log10(h_values)
log_error = np.log10(errors)

# Fit a straight line: log(error) = slope * log(h) + intercept
coeffs = np.polyfit(log_h, log_error, 1)
slope = coeffs[0]

# Plotting
plt.figure(figsize=(6, 4))
plt.plot(h_values, errors, 'o-', label=f'Slope ≈ {slope:.2f}')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Mesh size h (log scale)')
plt.ylabel('L2 error (log scale)')
plt.title('Log-Log Plot: L2 Error vs. Mesh Size')
plt.grid(True, which="both", ls="--")
plt.legend()
plt.tight_layout()
plt.show()
