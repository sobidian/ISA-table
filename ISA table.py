import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def atm_table():
    # Set up the altitude range
    alt = np.arange(0, 82000, 100)

    # Calculate the temperature
    temp = []
    for h in alt:
        if h <= 11000:
            T = 288.15 - 0.0065 * h
        elif h <= 25000:
            T = 216.65
        elif h <= 47000:
            T = 216.65 + 0.001 * (h - 25000)
        else:
            T = 282.65 + 0.0028 * (h - 47000)
        temp.append(T)

    # Calculate the pressure
    press = []
    for h in alt:
        if h <= 11000:
            p = 101325 * (288.15 / (288.15 - 0.0065 * h)) ** (9.80665 / (-0.0065 * 287))
        elif h <= 25000:
            p = 22632.1 * np.exp(-0.0001577 * (h - 11000))
        elif h <= 47000:
            p = 2480.62 * (216.65 / (216.65 + 0.001 * (h - 25000))) ** 34.1632
        else:
            p = 149.034 * ((282.65 + 0.0028 * (h - 47000)) / 282.65) ** (-34.1632 / 0.0065)
        press.append(p)

    # Calculate the density
    dens = []
    for i in range(len(alt)):
        rho = press[i] / (287 * temp[i])
        dens.append(rho)

    # Calculate relative temperature, pressure, and density
    T0 = temp[0]
    p0 = press[0]
    rho0 = dens[0]
    T_rel = np.array(temp) / T0
    p_rel = np.array(press) / p0
    rho_rel = np.array(dens) / rho0

    # Create a pandas dataframe for the table
    data = {'Altitude (m)': alt, 'Temperature (K)': temp, 'Pressure (Pa)': press, 'Density (kg/m^3)': dens,
            'Relative Temperature': T_rel, 'Relative Pressure': p_rel, 'Relative Density': rho_rel}
    table = pd.DataFrame(data)

    # Plot temperature, pressure, and density vs altitude
    fig, axs = plt.subplots(2, 3, figsize=(16, 8))

    axs[0, 0].plot(table['Temperature (K)'], table['Altitude (m)'])
    axs[0, 0].set_xlabel('Temperature (K)')
    axs[0, 0].set_ylabel('Altitude (m)')
    axs[0, 0].set_title('Temperature vs Altitude')

    axs[0, 1].plot(table['Pressure (Pa)'], table['Altitude (m)'])
    axs[0, 1].set_xlabel('Pressure (Pa)')
    axs[0, 1].set_ylabel('Altitude (m)')
    axs[0, 1].set_title('Pressure vs Altitude')

    axs[0, 2].plot(table['Density (kg/m^3)'], table['Altitude (m)'])
    axs[0, 2].set_xlabel('Density (kg/m^3)')
    axs[0, 2].set_ylabel('Altitude (m)')
    axs[0, 2].set_title('Density vs Altitude')

    axs[1, 0].plot(table['Relative Temperature'], table['Altitude (m)'])
    axs[1, 0].set_xlabel('Relative Temperature')
    axs[1, 0].set_ylabel('Altitude (m)')
    axs[1, 0].set_title('Relative Temperature vs Altitude')

    axs[1, 1].plot(table['Relative Pressure'], table['Altitude (m)'])
    axs[1, 1].set_xlabel('Relative Pressure')
    axs[1, 1].set_ylabel('Altitude (m)')
    axs[1, 1].set_title('Relative Pressure vs Altitude')

    axs[1, 2].plot(table['Relative Density'], table['Altitude (m)'])
    axs[1, 2].set_xlabel('Relative Density')
    axs[1, 2].set_ylabel('Altitude (m)')
    axs[1, 2].set_title('Relative Density vs Altitude')

    plt.tight_layout()
    plt.show()

    return table

# Call the function to generate the atmospheric properties table and plot
atm_table()


# In[ ]:




