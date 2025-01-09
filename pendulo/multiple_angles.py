import numpy as np
from pendulum_lib import pendulo_rk4, calcular_periodo, pendulo_rk4_adaptativo

def testar_varios_angulos():
    """
    Testa diferentes ângulos iniciais e passos constantes, exibindo os resultados.
    Também avalia o uso do passo adaptativo e compara com a solução analítica simplificada.
    """
    # Valores
    g = 9.81
    l = 10.0
    t_max = 10.0
    passos = [0.01, 0.001, 0.0001]
    angulos_iniciais = [np.radians(5), np.radians(20), np.radians(40), np.radians(60)]
    tolerancia_adaptativo = 10e-5

    resultados = []

    for theta0 in angulos_iniciais:
        for h in passos:
            t, theta = pendulo_rk4(g, l, theta0, 0.0, t_max, h)
            periodo_numerico = calcular_periodo(t, theta)
            periodo_analitico = 2 * np.pi * np.sqrt(l / g)
            resultados.append((
                np.degrees(theta0), h, "Fixo", periodo_numerico, periodo_analitico, abs(periodo_analitico - periodo_numerico)
            ))

        # Passo adaptativo
        t, theta = pendulo_rk4_adaptativo(g, l, theta0, 0.0, t_max, h_inicial=0.01, tolerancia=tolerancia_adaptativo)
        periodo_numerico = calcular_periodo(t, theta)
        periodo_analitico = 2 * np.pi * np.sqrt(l / g)
        resultados.append((
            np.degrees(theta0), "Adaptativo", tolerancia_adaptativo, periodo_numerico, periodo_analitico, abs(periodo_analitico - periodo_numerico)
        ))

    # Exibindo os resultados
    print("Ângulo Inicial (graus), Passo/Tolerância, Estratégia, Período Numérico, Período Analítico, Erro")
    for res in resultados:
        print(f"{res[0]:.2f}, {res[1]}, {res[2]}, {res[3]:.12f}, {res[4]:.12f}, {res[5]:.12f}")

if __name__ == "__main__":
    testar_varios_angulos()
