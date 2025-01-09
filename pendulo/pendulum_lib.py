import numpy as np
import matplotlib.pyplot as plt

def pendulo_rk4(g, l, theta0, omega0, t_max, h):
    """
    Simula o movimento de um pêndulo usando o método de Runge-Kutta de 4ª ordem.

    Parâmetros:
    - g: aceleração da gravidade (m/s²)
    - l: comprimento da haste (m)
    - theta0: ângulo inicial (rad)
    - omega0: velocidade angular inicial (rad/s)
    - t_max: tempo total de simulação (s)
    - h: passo de integração (s)

    Retorna:
    - t: vetor de tempos
    - theta: vetor de ângulos
    """
    # Função diferencial para o sistema
    def derivadas(t, y):
        theta, omega = y
        dtheta_dt = omega
        domega_dt = -(g / l) * np.sin(theta)
        return np.array([dtheta_dt, domega_dt])

    # Vetores para armazenar os resultados
    t = np.arange(0, t_max, h)
    theta = np.zeros(len(t))
    omega = np.zeros(len(t))

    # Condições iniciais
    theta[0] = theta0
    omega[0] = omega0

    # Método de Runge-Kutta de 4ª ordem
    for i in range(1, len(t)):
        y = np.array([theta[i-1], omega[i-1]])
        k1 = h * derivadas(t[i-1], y)
        k2 = h * derivadas(t[i-1] + h / 2, y + k1 / 2)
        k3 = h * derivadas(t[i-1] + h / 2, y + k2 / 2)
        k4 = h * derivadas(t[i-1] + h, y + k3)
        
        y_next = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        theta[i] = y_next[0]
        omega[i] = y_next[1]

    return t, theta

def calcular_periodo(t, theta):
    """
    Calcula o período do pêndulo monitorando mudanças de sinal na velocidade.

    Parâmetros:
    - t: vetor de tempos
    - theta: vetor de ângulos

    Retorna:
    - período médio
    """
    cruzamentos = []
    for i in range(1, len(theta)):
        if theta[i-1] * theta[i] < 0:
            cruzamentos.append(t[i-1] + (t[i] - t[i-1]) * abs(theta[i-1]) / (abs(theta[i-1]) + abs(theta[i])))

    periodos = [2 * (cruzamentos[i] - cruzamentos[i-1]) for i in range(1, len(cruzamentos))]
    return np.mean(periodos) if len(periodos) > 0 else None

def pendulo_rk4_adaptativo(g, l, theta0, omega0, t_max, h_inicial, tolerancia):
    """
    Simula o movimento de um pêndulo usando o método de Runge-Kutta de 4ª ordem
    com passo adaptativo.

    Parâmetros:
    - g: aceleração da gravidade (m/s²)
    - l: comprimento da haste (m)
    - theta0: ângulo inicial (graus)
    - omega0: velocidade angular inicial (rad/s)
    - t_max: tempo total de simulação (s)
    - h_inicial: passo de integração inicial (s)
    - tolerancia: tolerância para ajuste de passo

    Retorna:
    - t: vetor de tempos
    - theta: vetor de ângulos
    """
    # Função diferencial para o sistema
    def derivadas(t, y):
        theta, omega = y
        dtheta_dt = omega
        domega_dt = -(g / l) * np.sin(theta)
        return np.array([dtheta_dt, domega_dt])

    t = [0]
    theta = [theta0]
    omega = [omega0]
    h = h_inicial

    while t[-1] < t_max:
        y = np.array([theta[-1], omega[-1]])

        # Runge-Kutta de 4ª ordem com passo h
        k1 = h * derivadas(t[-1], y)
        k2 = h * derivadas(t[-1] + h / 2, y + k1 / 2)
        k3 = h * derivadas(t[-1] + h / 2, y + k2 / 2)
        k4 = h * derivadas(t[-1] + h, y + k3)
        
        y_next = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6

        # Estimativa para passo reduzido
        h_half = h / 2
        k1_half = h_half * derivadas(t[-1], y)
        k2_half = h_half * derivadas(t[-1] + h_half / 2, y + k1_half / 2)
        k3_half = h_half * derivadas(t[-1] + h_half / 2, y + k2_half / 2)
        k4_half = h_half * derivadas(t[-1] + h_half, y + k3_half)

        y_half_step = y + (k1_half + 2 * k2_half + 2 * k3_half + k4_half) / 6

        # Calcula o erro
        erro = np.linalg.norm(y_next - y_half_step)

        # Ajusta o passo
        if erro > tolerancia:
            h /= 2 
        elif erro < tolerancia / 2:
            h *= 2

        t.append(t[-1] + h)
        theta.append(y_next[0])
        omega.append(y_next[1])

    return np.array(t), np.array(theta)