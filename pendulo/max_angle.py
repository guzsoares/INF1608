import numpy as np
import matplotlib.pyplot as plt
from pendulum_lib import pendulo_rk4

def calcular_periodo_interpolacao(t, theta):
    """
    Calcula o período do pêndulo usando a interpolação linear para mudanças de sinal.

    Parâmetros:
    - t: vetor de tempos
    - theta: vetor de ângulos

    Retorna:
    - período médio
    """
    cruzamentos = []
    for i in range(1, len(theta)):
        if theta[i-1] * theta[i] < 0:
            v1, v2 = theta[i-1], theta[i]
            t1, t2 = t[i-1], t[i]
            cruzamento = t1 + (abs(v1) / (abs(v1) + abs(v2))) * (t2 - t1)
            cruzamentos.append(cruzamento)

    if len(cruzamentos) < 2:
        return None 

    periodos = [2 * (cruzamentos[i] - cruzamentos[i-1]) for i in range(1, len(cruzamentos))]
    return np.mean(periodos) if len(periodos) > 0 else None

def plotar_erro_angulo_maximo():
    """
    Plota o erro do período numérico em relação ao analítico para diferentes ângulos iniciais,
    destacando o ponto onde o erro excede o valor tolerado (0.001).
    """
    # Valores
    g = 9.81
    l = 10.0
    t_max = 10.0
    h = 0.001
    erro_tolerado = 0.001

    # Geração de ângulos iniciais em radianos
    angulos_iniciais = np.linspace(np.radians(1), np.radians(30), 300)
    erros = []
    angulos_graus = []

    for theta0 in angulos_iniciais:
        t, theta = pendulo_rk4(g, l, theta0, 0.0, t_max, h)
        periodo_numerico = calcular_periodo_interpolacao(t, theta)
        periodo_analitico = 2 * np.pi * np.sqrt(l / g)
        erro = abs(periodo_analitico - periodo_numerico)
        erros.append(erro)
        angulos_graus.append(np.degrees(theta0))

    # Plotando o gráfico
    plt.figure(figsize=(10, 6))
    plt.plot(angulos_graus, erros, label="Erro numérico vs analítico")
    plt.axhline(y=erro_tolerado, color='r', linestyle='--', label=f"Erro tolerado ({erro_tolerado})")
    plt.axvline(x=2.84, color='g', linestyle='--', label="Ângulo máximo permitido (~2.84°)")
    plt.xlabel("Ângulo Inicial (graus)")
    plt.ylabel("Erro Absoluto")
    plt.title("Erro do Período Numérico em Relação ao Analítico")
    plt.legend()
    plt.grid()
    plt.xlim(2, 5)  
    plt.ylim(0, 0.005)
    plt.show()

if __name__ == "__main__":
    plotar_erro_angulo_maximo()
