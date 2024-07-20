import csv
import numpy as np
import sys
from lmfit import Parameters, minimize, create_params
import matplotlib.pyplot as plt

num_act = 50
num_in = 100
num_ex = 50
start = 0
duration = 200
iterations = 12
days = 200

I0 = 100
A0 = 50
E0 = 50
S0 = 1000
R0 = 10

initial_conditions = [I0, E0, A0, S0, R0]

alpha = 1.5

beta = 0.001
c = 0.01

b = 0.1

mu = 11
gamma = 0.005
delta = 0.1
k = 0.1
tau = 0.001
epsilon = 0.00002
eta = 0.001

L1 = 2500
L2 = 17500
L3 = 140


def main():
    time = [0 for i in range(days)]
    active = [0 for i in range(days)]
    exhausted = [0 for i in range(days)]
    myeloma = [0 for i in range(days)]
    inactive = [0 for i in range(days)]
    cxcl9 = [0 for i in range(days)]
    resistant = [0 for i in range(days)]
    for i in range(iterations):
        with open(
                f'../output/act:{num_act}_in:{num_in}_ex:{num_ex}_start:{start}_duration:{duration}/CellCounts_{i}.csv',
                'r') as data:
            reader = csv.DictReader(data)
            for row in reader:
                # print(row["Time"])
                time[int(row["Time"])] += float(row["Time"])
                # print(time[int(row["Time"])])
                active[int(row["Time"])] += float(row["TCELLS"])
                exhausted[int(row["Time"])] += float(row[" EXTCELLS"])
                myeloma[int(row["Time"])] += float(row[" MYELOMA"])
                inactive[int(row["Time"])] += float(row[" Inactive Tcells"])
                resistant[int(row["Time"])] += float(row[" Resistant"])
                cxcl9[int(row["Time"])] += float(row[" CXCL9"])

    for i in range(days):
        time[i] /= iterations
        active[i] /= iterations
        inactive[i] /= iterations
        exhausted[i] /= iterations
        myeloma[i] /= iterations
        cxcl9[i] /= iterations
        resistant[i] /= iterations
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_yscale('log')

    plt.plot(time, active, 'g', label = 'Active T Cells')
    plt.plot(time, exhausted, 'r', label = 'Exhausted T Cells')
    plt.plot(time, myeloma, 'b', label = 'Myeloma')
    plt.plot(time, inactive, 'y', label = 'Inactive T Cells')
    plt.plot(time, cxcl9, 'k', label = 'CXCL9')
    plt.plot(time, resistant, 'r--', label = 'Resistant')
    plt.xlabel("Time (Days)")
    plt.ylabel("Number of Cells")
    plt.legend()
    plt.savefig(f'../output/figures/act:{num_act}_in:{num_in}_ex:{num_ex}_start:{start}_duration:{duration}')
    plt.show()

    # print(time)
    # print(active)
    # print(exhausted)
    # print(myeloma)
    # print(inactive)
    # print(cxcl9)

    # params = Parameters()
    # params.add('alpha', value=alpha, vary=True)#, min = 0, max = 0.2)
    # params.add('beta', value=beta, vary=False)
    # params.add('c', value=c, vary=False)
    # params.add('b', value=b, vary=False)
    # params.add('mu', value=mu, vary=False)
    # params.add('gamma', value=gamma, vary=False)
    # params.add('k', value=k, vary=False)
    # params.add('tau', value=tau, vary=False)
    # params.add('delta', value=delta, vary=False)
    # params.add('epsilon', value=epsilon, vary=False)
    # params.add('eta', value=eta, vary=False)
    # params.add('L1', value=L1, vary=False)
    # params.add('L2', value=L2, vary=False)
    # params.add('L3', value=L3, vary=False)

    # t = np.linspace(0, 200, 200 * 100)
    # t1 = np.linspace(0, 200, 200)

   
    # optimized = minimize(residuals, params,
    #                      method="leastsq", args=(t, inactive, exhausted, active, myeloma, resistant))
    # print("OPTIMIZED")
    # print(optimized.residual)

    # sol = euler_method(t, ODE, initial_conditions, optimized.params)
    # # sol = euler_method(t, ODE, initial_conditions, params)
    # print(sol)

    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1)
    # # ax.set_yscale('log')
    
    # plt.plot(t, sol[:, 0], 'b', label="I predicted")
    # plt.plot(t1, inactive, 'b--', label='I actual')
    
    # plt.plot(t, sol[:, 1], 'r', label="E predicted")
    # plt.plot(t1, exhausted, 'r--', label="E actual")
    
    # plt.plot(t, sol[:, 2], 'g', label="A predicted")
    # plt.plot(t1, active, 'g--', label="A actual")
    
    # plt.plot(t, sol[:, 3], 'y', label="S predicted")
    # plt.plot(t1, myeloma, 'y--', label="S actual")
    
    # plt.plot(t, sol[:, 4], 'k', label="R predicted")
    # plt.plot(t1, resistant, 'k--', label="R actual")
    
    # plt.legend()
    # plt.show()


def euler_method(t, ODE_system, initial_conditions, params):
    """
      Takes np array of linspace and does Euler's method
    """
    args = getParametersValues(params)

    sol = []

    step = (t[-1] - t[0]) / len(t)
    print(f"Step: {step}")

    y_n = initial_conditions
    # print(y_n[1])
    num_equations = len(y_n)

    for i in range(len(t)):

        sol.append(list(y_n))
        # print(y_n)
        # print(*ODE_system(*y_n, *args))
        for e in y_n:
            if np.isnan(e):
                sys.exit(1)

        for j in range(num_equations):

            # print(f"{j} Prime: {ODE_system(*y_n, *args)[j]}")
            # print(f"What E_Prime should be: {E_prime(y_n[1], y_n[3], y_n[4], y_n[2], beta, c)}")
            y_n[j] = ODE_system(*y_n, params)[j] * step + y_n[j]
            if abs(y_n[j]) < 0.1:
                y_n[j] = 0

    # print(sol)
    return np.array(sol)


def I_prime(I, alpha):
    return -I * alpha


def E_prime(E, S, R, A, beta, c):
    return beta * (S + R) * A - c * E


def A_prime(A, I, S, R, b, mu, gamma, alpha, L3):
    return b * I * alpha + mu * A * (1 - A / L3) - (A) * (S + R) * gamma


def S_prime(S, A, k, tau, L2):
    return k * S * (1 - S / L2) - A * tau * S


def R_prime(R, S, A, delta, epsilon, eta, L1):
    return delta * R * (1 - R / L1) - eta * A * R


def ODE(I, E, A, S, R, params):  # alpha, beta, c, b, mu, gamma, k, tau, delta, epsilon, eta, L1, L2, L3):

    return [
        I_prime(I, params['alpha']),
        E_prime(E, S, R, A, params["beta"], params["c"]),
        A_prime(A, I, S, R, params['b'], params['mu'], params['gamma'], params['alpha'], params['L3']),
        S_prime(S, A, params['k'], params['tau'], params['L2']),
        R_prime(R, S, A, params['delta'], params['epsilon'], params['eta'], params['L1'])]


def residuals(params , t, inactive, exhausted, active, sensitive, resistant):
    
    global initial_conditions
    sol = euler_method(t, ODE, initial_conditions, params)
    
    I_norm = []
    E_norm = []
    A_norm = []
    S_norm = []
    R_norm = []
    for i in range(200):
        I_norm.append(sol[i * 100, 0])
        E_norm.append(sol[i * 100, 1])
        A_norm.append(sol[i * 100, 2])
        S_norm.append(sol[i * 100, 3])
        R_norm.append(sol[i * 100, 4])
    I_res = np.absolute((np.array(I_norm) - np.array(inactive)).flatten())
    E_res = np.absolute((np.array(E_norm) - np.array(inactive)).flatten())
    A_res = np.absolute((np.array(A_norm) - np.array(inactive)).flatten())
    S_res = np.absolute((np.array(S_norm) - np.array(inactive)).flatten())
    R_res = np.absolute((np.array(R_norm) - np.array(inactive)).flatten())
    
    return (I_res + E_res + A_res + S_res + R_res).flatten()

    # return [vals['alpha'] ** 2 - 2 * vals['alpha'] + 2 + vals['beta'] ** 2 - 2 * vals['beta'] + 2, 0]


def getParametersValues(parameters_list):
    p = []
    for e in parameters_list:
        p.append(parameters_list[e].value)
    # print(p)
    return p


# E_prime(E, S, R, A, beta, c)




if (__name__ == "__main__"):
    main()
