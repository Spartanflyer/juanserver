from math import tan, cos, sin, acos, asin, atan, pi
import math


def main():
    # disp('agregar las siguientes variables en unidades de sistema internacional')
    h = []
    LL = []
    tetar = 0.66  # angulo de rodado
    tetas = 0.7854  # superior avalancha
    tetai = 0.68068  # inferior avalancha
    tetasd = 0.57596  # superior deslizamiento
    tetaid = 0.463786  # git initinferior deslizamiento
    afr = 0.5  # angulo de friccion
    densidad = 750  # densidad del material
    dp = 0.002  # diametro d eparticula
    n = 2  # revoluciones por segundo
    fi = 0.01745  # angulo deinclinaciondel material
    flujom = 0.033  # flujomasico
    L = 2  # longituddel cilindro
    dc = 0.15  # diametrodel cilindro
    # variablesinternas #

    R = dc / 2  # radiodel cilindro
    q = flujom / densidad
    r = dp / 2  # radiodeparticula

    # bucle for de 1 a longitud 2 del cilindro #
    N = 1000  # repeticiones
    deltat = L / N
    LL.insert(0, 0)

    Bd = (0.75 * flujom * tan(tetar)) / (pi * n * densidad * R ** 3)
    F0 = 1.75 * Bd ** 0.5  # formula paraangulodereposoentre 1 - 2
    F0 = 0.02095720906

    E0 = 0.06700853247

    L0 = 0.000087713943
    h.insert(0, R - (R ** 2 - (L0 / 2) ** 2) ** 0.5)

    # Cálculodealturamáximadedeslizamientointermitente
    alfa = acos((R - h[0]) / R)
    h1 = 0.002

    while (tan(tetasd) - (((3 * tan(afr)) * (alfa - (cos(alfa) * sin(alfa)))) / (2 * ((sin(alfa)) ** 3)))) > 0:
        h1 = h1 + 0.0000001
        alfa = acos((R - h1) / R)

    hd = h1  # Alturamáximadedeslizamiento

    for i in range(N):
        # calculo        centroides
        Xa = R * cos((asin((R - h[i]) / R)) - tetai)
        Ya = R * sin((asin((R - h[i]) / R)) - tetai)
        Xp = R * cos((asin((R - h[i]) / R)) - tetas)
        Yp = R * sin((asin((R - h[i]) / R)) - tetas)

        m1 = tan(tetai)  # Interseccióndeloscentroides
        m2 = tan(tetas)  # Interseccióndeloscentroides

        Xpsc = ((Yp - Ya) + (m1 * Xa) - (m2 * Xp)) / (m1 - m2)
        Ypsc = (m1 * (Xpsc - Xa)) + Ya

        Xabc = (Xa + Xp + Xpsc) / 3
        Yabc = (Ya + Yp + Ypsc) / 3

        s = (((Xabc - Xpsc) ** 2) + ((Yabc - Ypsc) ** 2)) ** 0.5
        neo = (atan((Yabc - Ypsc) / (Xabc - Xpsc)))

        # froudecalculado
        Fr = (n ** 2) * (R / 9.8)  # numero de froude
        gamma = tetas - tetai
        Frc = 0.5 * (R / s) * (900 * (gamma ** 2) / (pi) ** 2) * (
                    sin(neo) - tan(tetai) * cos(neo))  # numerodefroudedespejado

        if h[0] > hd:
            h.insert(i,
                     (1 / cos((tetasd + tetaid) / 2)) * (((3 * q * (tetasd - tetaid) * sin((tetasd + tetaid) / 2)) / (
                             ((2 * R - (h[i])) * (h[i])) ** 1.5 * 4 * pi * (
                                 2 * (1 - cos(tetasd - tetaid))) ** 0.5)) - fi))
            h.insert(i + 1, h[i] * deltat + h[i])
        elif Fr > Frc:
            h.insert(i, ((3 * q * sin(tetar)) / ((4 * pi * n * cos(tetar)) * (((2 * R) - (h[i])) * (h[i])) ** 1.5)) - (
                    fi / cos(tetar)))
            h.insert(i + 1, h[i] * deltat + h[i])

        else:
            h.insert(i, ((3 * q * sin(neo)) / ((4 * pi * n * cos(neo)) * (((2 * R) - (h[i])) * (h[i])) ** 1.5)) - (
                        fi / cos(neo)))
            h.insert(i + 1, h[i] * deltat + h[i])

    LL[i + 1] = i * deltat

    # plot(LL, h)


if __name__ == "__main__":
    main()
