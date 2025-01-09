import numpy as np
import time
from pendulum_lib import pendulo_rk4

def verificar_tempo_real():
    """
    Compara o tempo físico simulado com o tempo computacional da execução.
    """
    # Valores
    g = 9.81
    l = 10.0
    theta0 = np.radians(10)
    t_max = 10 * 2 * np.pi * np.sqrt(l / g)
    h = 0.01

    print(f"Simulando {t_max:.2f} segundos físicos com passo {h:.3f} segundos...")

    # Tempo real antes da execução
    start_time = time.time()

    # Simula o movimento
    t, theta = pendulo_rk4(g, l, theta0, 0.0, t_max, h)

    # Tempo real após a execução
    exec_time = time.time() - start_time

    # Resultados
    tempo_fisico = t[-1]
    print(f"Tempo físico simulado: {tempo_fisico:.2f} segundos")
    print(f"Tempo real de execução: {exec_time:.4f} segundos")

    if exec_time < tempo_fisico:
        print("Conclusão: O programa executa mais rápido que o tempo físico simulado.")
    else:
        print("Conclusão: O programa executa mais lento ou equivalente ao tempo físico simulado.")

if __name__ == "__main__":
    verificar_tempo_real()
