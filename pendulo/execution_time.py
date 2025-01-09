import numpy as np
import time
from pendulum_lib import pendulo_rk4, calcular_periodo, pendulo_rk4_adaptativo

def tempo_execucao_simulacao():
    """
    Mede o tempo de execução da simulação para 10 períodos,
    considerando diferentes estratégias de passo.
    """
    # Valores
    g = 9.81 
    l = 10.0
    theta0 = np.radians(10) 
    t_max = 10 * 2 * np.pi * np.sqrt(l / g)
    passos = [0.01, 0.001, 0.0001]
    tolerancia_adaptativo = 10e-5

    print("Estratégia, Passo/Tolerância, Tempo de Execução (s)")

    # Passos fixos
    for h in passos:
        start_time = time.time()
        pendulo_rk4(g, l, theta0, 0.0, t_max, h)
        exec_time = time.time() - start_time
        print(f"Passo Fixo, {h:.5f}, {exec_time:.4f}")

    # Passo adaptativo
    start_time = time.time()
    pendulo_rk4_adaptativo(g, l, theta0, 0.0, t_max, h_inicial=0.01, tolerancia=tolerancia_adaptativo)
    exec_time = time.time() - start_time
    print(f"Passo Adaptativo, Tolerância {tolerancia_adaptativo:.1e}, {exec_time:.4f}")

    # Solução analítica simplificada
    start_time = time.time()
    periodo_analitico = 2 * np.pi * np.sqrt(l / g)
    exec_time = time.time() - start_time
    print(f"Analítico, Simplificado, {exec_time:.4f}")

if __name__ == "__main__":
    tempo_execucao_simulacao()
