
# Giovanni Correra 09/2024 #

import numpy as np

# ------------------------------ Data ---------------------------------- #

T = -146.3 + 273.15  # Temperature in K
P = 3 + 1.01325  # Pressure in bar

x = np.array([0.9984, 0.0016, 0.00])  # Molar fractions of CH4, N2, CO2

Tc = np.array([190.4, 126.2, 304.1])  # Critical temperatures in K
Pc = np.array([46.0, 33.9, 73.8])  # Critical pressures in bar
om = np.array([0.011, 0.039, 0.239])  # Omega coefficient
MW = np.array([16.043, 28.013, 44.010])  # Molar masses in g/mol

# -------------------------- Data check -------------------------------- #

if not np.isclose(np.sum(x), 1):
    print('ERROR: feed molar fractions do not add up to 1')
    exit()


# -------------------------- Main script ------------------------------- #

def srk(T, P, Tc, Pc, om, y, phase):
    R = 8.3145
    RT = R * T
    RTc = R * Tc

    S = 0.48 + 1.574 * om - 0.176 * om ** 2
    k = (1 + S * (1 - np.sqrt(T / Tc))) ** 2
    a = (0.42748 * k * RTc ** 2) / Pc
    b = 0.08664 * RTc / Pc
    AS = a * P / (RT ** 2)
    BS = b * P / RT

    aM = np.sqrt(np.outer(a, a))
    bM = np.add.outer(b, b) / 2
    am = np.dot(y, np.dot(aM, y))
    bm = np.dot(y, np.dot(bM, y))

    A = am * P / (RT ** 2)
    B = bm * P / RT

    alfa = -1
    beta = A - B - B ** 2
    gamma = -A * B

    # Analytical solution
    p = beta - (alfa ** 2) / 3
    q = 2 * (alfa ** 3) / 27 - alfa * beta / 3 + gamma
    q2 = q / 2
    a3 = alfa / 3
    D = (q ** 2) / 4 + p ** 3 / 27

    if D > 0:
        Z1 = np.cbrt(-q2 + np.sqrt(D)) + np.cbrt(-q2 - np.sqrt(D)) - a3
        Z = [Z1, Z1, Z1]
    elif D == 0:
        Z1 = -2 * np.cbrt(q2) - a3
        Z2 = np.cbrt(q2) - a3
        Z = [Z1, Z2, Z2]
    else:
        r = np.sqrt(-p ** 3 / 27)
        theta = np.arccos(-q2 * np.sqrt(-27 / p ** 3))
        Z1 = 2 * np.cbrt(r) * np.cos(theta / 3) - a3
        Z2 = 2 * np.cbrt(r) * np.cos((2 * np.pi + theta) / 3) - a3
        Z3 = 2 * np.cbrt(r) * np.cos((4 * np.pi + theta) / 3) - a3
        Z = [Z1, Z2, Z3]

    Z = np.sort(Z)

    if phase == 1:
        Z = np.max(Z)
    elif phase == 2:
        Z = np.min(Z)

    return Z


def density(Z, T, P, x, MW):
    v = Z * T * 8.3451 / (P * 1e5)  # Molar volume in m3/mol
    MW_tot = np.sum(x * MW)  # Molecular weight of mixture in g/mol
    rho = (1 / v) * MW_tot / 1e3  # Mixture density in kg/m3
    return rho


Z = srk(T, P, Tc, Pc, om, x, phase=2)
rho = density(Z, T, P, x, MW)

# ------------------------ Post - Processing --------------------------- #

print(f'Liquid density = {rho:.3f} [kg/m3]')
print()
print(f'T = {T - 273.15:.3f} [°C]    ', end='')
print(f'P, abs = {P:.5f} [bar]    ', end='')
print(f'P, rel = {P - 1.01325:.5f} [bar]')
print()
print('Liquid composition [mol/mol]')
print(f'x, CH4 = {x[0] * 100:.4f} %')
print(f'x, N2  = {x[1] * 100:.4f}  %')
print(f'x, CO2 = {x[2] * 100:.4f}  %')
print()
