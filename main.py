import numpy as np
from set_stack import set_stack
from ATR1D import ATR1D
from scipy.optimize import differential_evolution
import random
import time

# ... (keep the existing input data)

def calculate_TR_avg(thicknesses):
    start_time = time.time()
    try:
        dPSC = list(thicknesses[:8]) + [dgls] + list(thicknesses[8:])
        materials = ["air"] + matPSC + ["SiO2.csv", "EVA.csv", "const=1.24", "air"]
        d = np.array(dPSC + [50, 100, 200])
        
        stack = set_stack(materials, d, lam, theta, incoh)
        R_PSC_AZO_EVA_PV, T_PSC_AZO_EVA_PV, _ = ATR1D(stack)
        T_avg = np.mean(T_PSC_AZO_EVA_PV['sp'])
        R_avg = np.mean(R_PSC_AZO_EVA_PV['sp'])
        end_time = time.time()
        print(f"calculate_TR_avg took {end_time - start_time:.2f} seconds")
        return T_avg, R_avg
    except Exception as e:
        print(f"Error in calculate_TR_avg: {str(e)}")
        print(f"Input thicknesses: {thicknesses}")
        raise

def monte_carlo_optimization(num_iterations=100):
    start_time = time.time()
    best_thicknesses = dBragg1 + dBragg2
    best_T_avg, best_R_avg = calculate_TR_avg(best_thicknesses)
    progress = []

    print(f"Starting Monte Carlo optimization with {num_iterations} iterations")
    print(f"Initial T_avg: {best_T_avg:.6f}, R_avg: {best_R_avg:.6f}")

    for i in range(num_iterations):
        iteration_start = time.time()
        thicknesses = [random.uniform(10, 200) for _ in range(16)]
        T_avg, R_avg = calculate_TR_avg(thicknesses)

        if T_avg > best_T_avg:
            best_thicknesses = thicknesses
            best_T_avg = T_avg
            best_R_avg = R_avg
            print(f"Iteration {i+1}: New best T_avg = {best_T_avg:.6f}, R_avg = {best_R_avg:.6f}")

        progress.append((best_T_avg, best_R_avg))
        iteration_end = time.time()
        print(f"Iteration {i+1} took {iteration_end - iteration_start:.2f} seconds")

        if (i+1) % 10 == 0:
            print(f"Completed {i+1} iterations. Current best T_avg = {best_T_avg:.6f}, R_avg = {best_R_avg:.6f}")

    end_time = time.time()
    print(f"Monte Carlo optimization took {end_time - start_time:.2f} seconds")
    return best_thicknesses, best_T_avg, best_R_avg, progress

def genetic_algorithm_optimization(pop_size=20, num_generations=50):
    start_time = time.time()
    bounds = [(10, 200) for _ in range(16)]  # 16 layers to optimize
    
    progress = []
    def callback(xk, convergence):
        T_avg, R_avg = calculate_TR_avg(xk)
        progress.append((T_avg, R_avg))
        if len(progress) % 10 == 0:
            print(f"Generation {len(progress)}: T_avg = {T_avg:.6f}, R_avg = {R_avg:.6f}")

    print(f"Starting Genetic Algorithm optimization with population size {pop_size} and {num_generations} generations")
    
    try:
        result = differential_evolution(
            objective_function,
            bounds,
            popsize=pop_size,
            maxiter=num_generations,
            mutation=(0.5, 1.5),
            recombination=0.7,
            updating='deferred',
            workers=-1,  # Use all available CPU cores
            callback=callback,
            tol=1e-10,
            polish=True
        )
        final_T_avg, final_R_avg = calculate_TR_avg(result.x)
        end_time = time.time()
        print(f"Genetic Algorithm optimization took {end_time - start_time:.2f} seconds")
        return result.x, final_T_avg, final_R_avg, progress
    except Exception as e:
        print(f"Error in genetic_algorithm_optimization: {str(e)}")
        raise

print("Calculating initial T_avg and R_avg...")
initial_T_avg, initial_R_avg = calculate_TR_avg(dBragg1 + dBragg2)
print(f"Initial T_avg: {initial_T_avg:.4f}, R_avg: {initial_R_avg:.4f}")

print("\nRunning Monte Carlo optimization...")
mc_thicknesses, mc_T_avg, mc_R_avg, mc_progress = monte_carlo_optimization()
print(f"Monte Carlo optimization result: T_avg = {mc_T_avg:.6f}, R_avg = {mc_R_avg:.6f}")
print(f"Optimized thicknesses: {[f'{t:.2f}' for t in mc_thicknesses]}")

print("\nRunning Genetic Algorithm optimization...")
ga_thicknesses, ga_T_avg, ga_R_avg, ga_progress = genetic_algorithm_optimization()
print(f"Genetic Algorithm optimization result: T_avg = {ga_T_avg:.6f}, R_avg = {ga_R_avg:.6f}")
print(f"Optimized thicknesses: {[f'{t:.2f}' for t in ga_thicknesses]}")

print("\nOptimization complete.")
