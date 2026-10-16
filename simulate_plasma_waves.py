from scipy.fft import fft, ifft, fftfreq
from numpy import pi, linspace, sin, cos, zeros_like, concatenate
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# *_hat will be used to denote normalised variables
# *_bar will be used to denote their Fourier transforms

def calculate_electric_field(n_e_hat, x_hat):
    # Calculate total charge density using electron number density
    # rho = q_p(n_p_bar - n_e_bar)
    # Normalised such that q_p and n_p_bar = 1
    # rho = (1 - n_e_bar)
    rho_hat = 1. - n_e_hat;

    # Calculate the Fourier transform of total charge density to get it in terms of frequency
    rho_bar = fft(rho_hat);

    # Get the samples per unit spacing (wavenumber) of the ft
    N = len(x_hat);
    dx = (x_hat[-1] - x_hat[0]) / (N - 1);
    wavenums = fftfreq(N, dx);

    # Only the first half of these results contains the positive wavenumbers we want
    # Defining the last index of interest
    upperlim = N//2;

    # Calculate the fourier transform of E_bar from the ft of rho_bar
    # Assume that the average electric field is zero - so leaves E_bar[0] = 0
    E_bar = zeros_like(wavenums, dtype=complex);
    E_bar[1:] = rho_bar[1:] / (1j * wavenums[1:]);

    # Reverse ft on E_bar to get to E_hat (spatial domain)
    return ifft(E_bar).real;

def fourier_derivative ( field , xx ):
    """ This function calculates d field / dx using a Fourier approach """
    wavenumbers = fftfreq ( len ( xx ) , xx [1] - xx [0])
    return ifft ( complex (0 , 1) * wavenumbers * fft ( field )). real

# Define function d[n_e_hat, U_e_hat]/dt_hat = [-(U_e_hat * dn_e_hat/dx_hat + n_e_hat * dU_e_hat/dx_hat), (E_x_hat*q_e_hat)/m_e_hat - T_e_hat/(m_e_hat*n_e_hat) * dn_e_hat/dx_hat]
def dfdt(curruent_t, current_f, x_hat, T_e_hat):
    # Define electron constants
    q_e_hat = -1.;
    # m_e_hat = 1;

    N = len(x_hat);
    n_e_hat = current_f[:N];
    U_e_hat = current_f[N:];

    # Calculate spacial derivatives
    dn_e_hat_dx_hat = fourier_derivative(n_e_hat, x_hat);
    dU_e_hat_dx_hat = fourier_derivative(U_e_hat, x_hat);

    # Calculate the E_hat for the current n_e_hat
    E_x_hat = calculate_electric_field(n_e_hat, x_hat);

    # Calculate temporal derivatives
    # m_e_hat has been normalised to 1
    dn_e_hat_dt_hat = -(U_e_hat * dn_e_hat_dx_hat + n_e_hat * dU_e_hat_dx_hat);
    dU_e_hat_dt_hat = (E_x_hat * q_e_hat) - (T_e_hat/n_e_hat) * dn_e_hat_dx_hat;

    return concatenate([dn_e_hat_dt_hat, dU_e_hat_dt_hat]);

def initial_ne(xx, kx, amp):
    nx = len(xx);
    ne0 = 1 + amp * cos(2*pi * xx * kx);
    ne0[:nx//4] = 1.0;
    ne0[3*nx//4:] = 1.0;
    return ne0

def initial_Ue ( xx ):
    return xx * 0.0

def simulate_electron_density_in_spacetime(x_hat, t_hat, T_e_hat, method="RK45"):
    # Create initial values
    N = len(x_hat);
    n_e_hat0 = initial_ne(x_hat, 4., 0.01);
    U_e_hat0 = initial_Ue(x_hat);

    y0 = concatenate([n_e_hat0, U_e_hat0]);

    result = solve_ivp(
        fun=dfdt,
        t_span=[t_hat[0], t_hat[-1]],
        y0=y0,
        t_eval=t_hat,
        args=[x_hat, T_e_hat],
        method=method
    );

    return result.y[:N], result.y[N:];

def plot_electron_density_in_spacetime(t_hat, x_hat, n_e_hats, title):
    pcm = plt.pcolormesh(t_hat, x_hat, n_e_hats);
    plt.colorbar(pcm, label=r"Electron density normalised as $\hat{n}_e = \frac{n_e}{n_p}$");
    plt.xlabel(r"Time normalised as $\hat{t} = t \omega_p$");
    plt.ylabel(r"Distance normalised as $\hat{x} = \frac{x}{\lambda_D}$");
    plt.title(title);
    plt.show();

def main():
    # Calculate and plot electric field
    # ---------------------------------
    x_hat = linspace(0, 1, 101, endpoint=False);
    n_e_hat = 1 + sin(2*pi * x_hat);

    E_hat = calculate_electric_field(n_e_hat, x_hat);

    one_minus_dE_hat_dx_hat = 1 - fourier_derivative(E_hat, x_hat);

    plt.plot(x_hat, n_e_hat, color="orange", label=r"Electron density normalised as $\hat{n}_e = \frac{n_e}{n_p}$");
    plt.plot(x_hat, E_hat,color="red", label=r"Electric field normalised as $\hat{E}_x = \frac{E_x \lambda_D}{T_p}$");
    plt.plot(x_hat, one_minus_dE_hat_dx_hat, linestyle="dotted", color="green", label=r"Recalculated normalised electron density");
    plt.xlabel(r"Distance normalised as $\hat{x} = \frac{x}{\lambda_D}$");
    plt.ylabel("Normalised values of electron density and electric field");
    plt.title(
        r"""
        Normalised electric field, $\hat{E}_x$, of a plasma with normalised electron density, $\hat{n}_e$.
        $\hat{E}_x$ verified by comparing initial $\hat{n}_e$ against $\hat{n}_e = 1 - \frac{d \hat{E}_x}{d \hat{x}}$
        """
        );
    plt.grid();
    plt.legend();
    plt.show();

    print("Please wait - Solving IVPs...");
    # Finding the time derivatives T_e_hat = 0
    # ----------------------------------------
    t_hat_T0 = linspace(0, 2*pi, 251, endpoint=False);
    x_hat_T0 = linspace(0, 1, 401, endpoint=False);
    T_e_hat_T0=0.
    title_T0 = r"Evolution of normalised electron density, $\hat{n}_e$, in normalised time, $\hat{t}$." + "\n" + r"Electron temperature normalised as $\hat{T}_e = \frac{T_e}{T_p}$ is equal to " + str(T_e_hat_T0) + ". Using RK45 method to solve IVP.";

    n_e_hats_T0, _ = simulate_electron_density_in_spacetime(x_hat_T0, t_hat_T0, T_e_hat_T0);

    # Finding the time derivatives T_e_hat = 1 using RK45 method
    # ----------------------------------------------------------
    t_hat_RK45 = linspace(0, 2*pi, 251, endpoint=False);
    x_hat_RK45 = linspace(0, 1, 401, endpoint=False);
    T_e_hat_RK45=1.
    title_RK45 = r"Evolution of normalised electron density, $\hat{n}_e$, in normalised time, $\hat{t}$." + "\n" + r"Electron temperature normalised as $\hat{T}_e = \frac{T_e}{T_p}$ is equal to " + str(T_e_hat_RK45) + ". Using RK45 method to solve IVP.";

    n_e_hats_RK45, _ = simulate_electron_density_in_spacetime(x_hat_RK45, t_hat_RK45, T_e_hat_RK45);

    # Finding the time derivatives T_e_hat = 1 using Radau method
    # -----------------------------------------------------------
    t_hat_Rad = linspace(0, 2*pi, 251, endpoint=False);
    x_hat_Rad = linspace(0, 1, 401, endpoint=False);
    T_e_hat_Rad=1.
    title_Rad = r"Evolution of normalised electron density, $\hat{n}_e$, in normalised time, $\hat{t}$." + "\n" + r"Electron temperature normalised as $\hat{T}_e = \frac{T_e}{T_p}$ is equal to " + str(T_e_hat_Rad) + ". Using Radau method to solve IVP.";

    n_e_hats_Rad, _ = simulate_electron_density_in_spacetime(x_hat_Rad, t_hat_Rad, T_e_hat_Rad, method="Radau");

    # Produce all the plots
    # ---------------------
    plot_electron_density_in_spacetime(t_hat_T0, x_hat_T0, n_e_hats_T0, title_T0);
    plot_electron_density_in_spacetime(t_hat_RK45, x_hat_RK45, n_e_hats_RK45, title_RK45);
    plot_electron_density_in_spacetime(t_hat_Rad, x_hat_Rad, n_e_hats_Rad, title_Rad);

    plt.plot(x_hat_RK45, n_e_hats_RK45[:, 199], label="RK45");
    plt.plot(x_hat_Rad, n_e_hats_Rad[:, 199], label="Radau");
    plt.legend();
    plt.title("Comparison of the accuracy of RK45 and Radau methods for solving electron density evolution in time.\n" + r"Slice in time taken at $\hat{t}$ = " + str("%.2g" % t_hat_Rad[199]));
    plt.xlabel(r"Distance normalised as $\hat{x} = \frac{x}{\lambda_D}$");
    plt.ylabel(r"Electron density normalised as $\hat{n}_e = \frac{n_e}{n_p}$");
    plt.show();

if __name__ == "__main__":
    main();
