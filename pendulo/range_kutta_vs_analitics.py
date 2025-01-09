import numpy as np
import matplotlib.pyplot as plt
from pendulum_lib import pendulo_rk4_adaptativo

def comparar_rk_adaptativo_com_analitico():
    """
    Compara o método de Runge-Kutta adaptativo com a solução analítica simplificada
    para ângulos iniciais específicos.
    """
    # Valores
    g = 9.81
    l = 10.0 
    t_max = 10.0
    h_inicial = 0.01
    tolerancia = 10e-5
    angulos_iniciais = [5, 20, 40, 60]

    print("Ângulo, Período Numérico, Período Analítico, Erro Relativo")

    for angulo in angulos_iniciais:
        theta0 = np.radians(angulo)

        # Solução numérica com RK adaptativo
        t, theta = pendulo_rk4_adaptativo(g, l, theta0, 0.0, t_max, h_inicial, tolerancia)
        cruzamentos = []
        for i in range(1, len(theta)):
            if theta[i-1] * theta[i] < 0:
                cruzamento = t[i-1] + (t[i] - t[i-1]) * abs(theta[i-1]) / (abs(theta[i-1]) + abs(theta[i]))
                cruzamentos.append(cruzamento)

        if len(cruzamentos) > 1:
            periodos = [2 * (cruzamentos[i] - cruzamentos[i-1]) for i in range(1, len(cruzamentos))]
            periodo_rk_adaptativo = np.mean(periodos)
        else:
            periodo_rk_adaptativo = None

        # Solução analítica simplificada
        periodo_analitico = 2 * np.pi * np.sqrt(l / g)

        # Erro relativo
        if periodo_rk_adaptativo is not None:
            erro_relativo = abs(periodo_analitico - periodo_rk_adaptativo) / periodo_analitico
        else:
            erro_relativo = None

        # Exibir resultados
        print(f"{angulo:.2f}, {periodo_rk_adaptativo:.10f}, {periodo_analitico:.10f}, {erro_relativo:.10f}")

def gerar_graficos(angulos_iniciais = [np.radians(5), np.radians(20), np.radians(40), np.radians(60)]):
    """
    Gera gráficos comparativos entre as soluções analítica e numérica
    para diferentes ângulos iniciais.
    """
    # Valores
    g = 9.81 
    l = 10.0
    t_max = 10.0
    h = 0.01

    for theta0 in angulos_iniciais:
        t, theta = pendulo_rk4_adaptativo(g, l, theta0, 0.0, t_max, h_inicial=0.01, tolerancia=10e-5)
        theta_analitico = theta0 * np.cos(np.sqrt(g / l) * t)

        plt.figure()
        plt.plot(t, theta, label="RK4 (h=0.01)")
        plt.plot(t, theta_analitico, label="Solução Analítica", linestyle="dashed")
        plt.xlabel("Tempo (s)")
        plt.ylabel("Ângulo (rad)")
        plt.legend()
        plt.title(f"Movimento do Pêndulo - Ângulo Inicial: {np.degrees(theta0):.2f}°")
        plt.grid()
        plt.show()

if __name__ == "__main__":
    comparar_rk_adaptativo_com_analitico()
    gerar_graficos()
