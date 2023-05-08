import random

def monte_carlo_ln(t, n):
    """Approximates the natural logarithm of t using Monte Carlo integration
    with n random points."""
    total = 0
    for i in range(n):
        x = random.uniform(1, t)
        total += 1 / x
    return (t - 1) * total / n

# Example usage:
t = 5
n = 1000000
approx_ln = monte_carlo_ln(t, n)
print(f"ln({t}) ≈ {approx_ln:.6f}")
