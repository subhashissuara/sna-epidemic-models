# --------------------------
#   Author: Subhashis Suara
#   Roll No: UCSE19012
# --------------------------

import numpy as np
import networkx as nx
import matplotlib.pylab as plt
from scipy.integrate import odeint

# ------------- EDIT BELOW THIS LINE -------------
N = 1000
k = 4
p = 0.005
S = N - 1
I = 1
R = 0
beta = 0.6
gamma = 0.05
mu = 0.2
lambdaa = 0.4
number_of_days = 180
# ------------- EDIT ABOVE THIS LINE -------------

def create_graph():
    global N, k, p
    G = nx.watts_strogatz_graph(N, k, p)
    return G

def si_model(G):
    global N, S, I, beta

    S_local = S
    I_local = I
    susceptible = []
    infected = []

    for time in range (0, number_of_days):
        S_local, I_local = S_local - beta * ((S_local * I_local / N)), I_local + beta * ((S_local * I_local) / N)
        susceptible.append(S_local)
        infected.append(I_local)

    figure = plt.figure()
    figure.canvas.set_window_title('SI model')

    susceptible_line, = plt.plot(susceptible, label='Susceptible(t)')
    infected_line, = plt.plot(infected, label='Infected(t)')

    plt.legend(handles=[susceptible_line, infected_line])
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('Time in Days')
    plt.ylabel('Population Fraction')
    plt.show()

def sis_model(G):
    global N, S, I, beta, gamma

    S_local = S
    I_local = I
    susceptible = []
    infected = []

    for time in range (0, number_of_days):
        S_local, I_local = S_local - (((beta * S_local * I_local) / N) + (gamma * I_local)), I_local + (((beta * S_local * I_local) / N) - (gamma * I))
        susceptible.append(S_local)
        infected.append(I_local)

    figure = plt.figure()
    figure.canvas.set_window_title('SIS model')

    susceptible_line, = plt.plot(susceptible, label='Susceptible(t)')
    infected_line, = plt.plot(infected, label='Infected(t)')

    plt.legend(handles=[susceptible_line, infected_line])
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('Time in Days')
    plt.ylabel('Population Fraction')
    plt.show()

def sir_diff(sir, t):
    global N, S, I, R, beta, gamma

    dsdt = - (beta * sir[0] * sir[1]) / N
    didt = (beta * sir[0] * sir[1]) / N - gamma * sir[1]
    drdt = gamma * sir[1]
    dsirdt = [dsdt, didt, drdt]
    return dsirdt

def sir_model(G):
    global N, S, I, R, beta, gamma

    sir_initial = (S, I, R)
    time = np.linspace(0, number_of_days)
    sir = odeint(sir_diff, sir_initial, time)

    figure = plt.figure()
    figure.canvas.set_window_title('SIR model')

    plt.plot(time, sir[:, 0], label='Susceptible(t)')
    plt.plot(time, sir[:, 1], label='Infected(t)')
    plt.plot(time, sir[:, 2], label='Recovered(t)')
    plt.legend()
    plt.xlabel('Time in Days')
    plt.ylabel('Population Fraction')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.show()

def sirs_model(G):
    global N, S, I, R, beta, gamma, lambdaa, mu

    susceptible = []
    infected = []
    recovered = []
    
    for time in range (1, number_of_days):
        S, I = S - (beta * S * I) / N + lambdaa * R, I + ((beta * S * I) / N) - lambdaa * R
        R = mu * I - lambdaa * S
        susceptible.append(S)
        infected.append(I)
        recovered.append(R)

    figure = plt.figure()
    figure.canvas.set_window_title('SIRS model')

    susceptible_line, = plt.plot(susceptible, label='Susceptible(t)')
    infected_line, = plt.plot(infected, label='Infected(t)')
    recovered_line, = plt.plot(recovered, label='Recovered(t)')

    plt.legend(handles=[susceptible_line, infected_line, recovered_line])
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('Time in Days')
    plt.ylabel('Population Fraction')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.show()

def main():
    G = create_graph()
    si_model(G)
    sis_model(G)
    sir_model(G)
    sirs_model(G)

if __name__ == '__main__':
    main()